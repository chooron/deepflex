"""
$(TYPEDEF)
A struct representing common hydrological fluxes
# Fields
$(FIELDS)
# Example
```
# define a common function
function flow_func(
    i::namedtuple(:baseflow, :surfaceflow),
    p::NamedTuple;
    kw...
)
    i[:baseflow] .+ i[:surfaceflow]
end

flow_flux = SimpleFlux(
    [:baseflow, :surfaceflow],
    :flow,
    param_names=Symbol[],
    func=flow_func
)
```
"""
struct SimpleFlux <: AbstractSimpleFlux
    "input variables information"
    input_info::NamedTuple
    "output variables information"
    output_info::NamedTuple
    "parameters information"
    param_info::NamedTuple
    """
    functions used to implement hydrological flux computation,
    It requires the input format to be (i::Vector, p::Vector)
    """
    inner_func::Function
    "flux expressions"
    flux_exprs::Vector{Num}

    function SimpleFlux(
        input_info::NamedTuple,
        output_info::NamedTuple,
        param_info::NamedTuple,
        inner_func::Function,
        flux_exprs::Vector{Num}
    )
        return new(
            input_info,
            output_info,
            param_info,
            inner_func,
            flux_exprs
        )
    end

    function SimpleFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        param_names::Vector{Symbol}=Symbol[];
        flux_funcs::Vector{<:Function}=Function[],
        kwargs...
    )
        #* 获取输入输出名称
        input_names, output_names = flux_names[1], flux_names[2]

        if length(flux_funcs) > 0
            inputs = [first(@variables $var(t) = 0.0) for var in input_names]
            outputs = [first(@variables $var(t) = 0.0) for var in output_names]
            params = [first(@parameters $var = 0.0) for var in param_names]
            flux_exprs = [flux_func(inputs, params) for flux_func in flux_funcs]
        else
            #* 根据输入输出参数名称获取对应的计算公式
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            flux_flux_exprs = expr(hydro_equation; kwargs...)
            #* 得到计算函数
            flux_funcs = [
                build_function(hydro_expr, hydro_equation.inputs, hydro_equation.params, expression=Val{false})
                for hydro_expr in flux_flux_exprs
            ]
            #* 得到计算公式
            flux_exprs = [hydro_expr for hydro_expr in flux_flux_exprs]
            inputs, outputs, params = hydro_equation.inputs, hydro_equation.outputs, hydro_equation.params
        end

        function inner_flux_func(input::AbstractVector, params::AbstractVector)
            [func(input, params) for func in flux_funcs]
        end

        function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
            reduce(hcat, [func.(eachrow(input), Ref(params)) for func in flux_funcs])
        end

        return new(
            NamedTuple{Tuple(input_names)}(inputs),
            NamedTuple{Tuple(output_names)}(outputs),
            NamedTuple{Tuple(param_names)}(params),
            inner_flux_func,
            flux_exprs,
        )
    end

    function SimpleFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        flux_exprs::Vector{Num}
    )
        inputs, outputs = fluxes[1], fluxes[2]
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params)

        #* 得到计算函数
        flux_funcs = [build_function(flux_expr, inputs, params, expression=Val{false}) for flux_expr in flux_exprs]

        function inner_flux_func(input::AbstractVector, params::AbstractVector)
            [func(input, params) for func in flux_funcs]
        end

        function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
            reduce(hcat, [func.(eachrow(input), Ref(params)) for func in flux_funcs])
        end

        return new(
            NamedTuple{(Tuple(input_names))}(inputs),
            NamedTuple{(Tuple(output_names))}(outputs),
            NamedTuple{(Tuple(param_names))}(params),
            inner_flux_func,
            flux_exprs,
        )
    end
end

struct StateFlux <: AbstractStateFlux
    "input variables information"
    input_info::NamedTuple
    "output variables information"
    output_info::NamedTuple
    "parameters information"
    param_info::NamedTuple
    "state expression"
    state_expr::Num

    function StateFlux(
        fluxes::Vector{Num},
        state::Num,
        params::Vector{Num};
        state_expr::Num
    )
        #* 转换为Symbol
        state_input_names = Symbolics.tosymbol.(fluxes, escape=false)
        state_name = Symbolics.tosymbol(state, escape=false)
        state_param_names = Symbolics.tosymbol(params, escape=false)
        state_input_ntp = NamedTuple{Tuple(state_input_names)}(fluxes)
        state_param_ntp = NamedTuple{Tuple(state_param_names)}(params)

        return new(
            state_input_ntp,
            NamedTuple{tuple(state_name)}([state]),
            state_param_ntp,
            state_expr,
        )
    end

    function StateFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        state::Num;
    )
        influxes, outfluxes = fluxes[1], fluxes[2]
        #* 构建函数和公式
        state_expr = sum(influxes) - sum(outfluxes)
        return StateFlux(vcat(influxes, outfluxes), state, Num[], state_expr=state_expr)
    end

    function StateFlux(
        states::Pair{Num,Num};
    )
        #* 构建函数和公式
        ori_state, new_state = states[1], states[2]
        state_expr = new_state - ori_state
        return StateFlux([new_state], ori_state, Num[], state_expr=state_expr)
    end

    function StateFlux(
        flux_names::Vector{Symbol},
        state_name::Symbol,
        param_names::Vector{Symbol};
        state_func::Function
    )
        #* 转换为Symbol
        fluxes = [first(@variables $nm(t) = 0.0) for nm in flux_names]
        state = first(@variables $state_name(t) = 0.0)
        params = [first(@parameters $nm = 0.0) for nm in param_names]
        state_expr = state_func(fluxes, params)
        return StateFlux(fluxes, state, params, state_expr=state_expr)
    end

    function StateFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        state_name::Symbol;
    )
        influx_names, outflux_names = flux_names[1], flux_names[2]
        state_input_names = vcat(influx_names, outflux_names)
        #* 构建函数和公式
        state_func = (i, p) -> sum(i[1:length(influx_names)]) - sum(i[length(influx_names)+1:length(influx_names)+length(outflux_names)])
        return StateFlux(state_input_names, state_name, Symbol[], state_func=state_func)
    end

    function StateFlux(
        state_names::Pair{Symbol,Symbol};
    )
        ori_state_name, new_state_name = state_names[1], state_names[2]
        #* 构建函数和公式
        ori_state = first(@variables $ori_state_name(t) = 0.0)
        new_state = first(@variables $new_state_name(t) = 0.0)
        return StateFlux(ori_state => new_state)
    end
end


"""
$(TYPEDEF)
A flux used in hydrological unit-hydrograph calculations
# Fields
$(FIELDS)
# Example
```
slowflow_lagflux = LagFlux(:slowflow, :slowflow_lag, lag_func=LumpedHydro.uh_1_half, lag_time=:x4)
```
"""
struct LagFlux <: AbstractLagFlux
    "input variables information"
    input_info::NamedTuple
    "output variables information"
    output_info::NamedTuple
    "parameters information"
    param_info::NamedTuple
    "function used to represent hydrograph unit"
    inner_func::Function

    function LagFlux(
        flux_names::Pair{Symbol,Symbol},
        lag_time_name::Symbol,
        lag_func::Function;
        kwargs...,
    )
        function inner_flux_func(i::AbstractArray, p::AbstractVector)
            lag_flux = solve_lag_flux(i, p[1], lag_func, kwargs...)
            [lag_flux]
        end

        fluxes = [first(@variables $nm(t) = 0.0) for nm in flux_names]
        lag_time = first(@variables $lag_time_name = 0.0)

        new(
            NamedTuple{(flux_names[1],)}([fluxes[1]]),
            NamedTuple{(flux_names[2],)}([fluxes[2]]),
            NamedTuple{(lag_time_name,)}([lag_time]),
            inner_flux_func,
        )
    end

    function LagFlux(
        fluxes::Pair{Num,Num},
        lag_time::Num,
        lag_func::Function;
        kwargs...,
    )
        function inner_flux_func(i::AbstractArray, p::AbstractVector)
            lag_flux = solve_lag_flux(i, p[1], lag_func, kwargs...)
            [lag_flux]
        end

        flux_name_1 = Symbolics.tosymbol(fluxes[1], escape=false)
        flux_name_2 = Symbolics.tosymbol(fluxes[2], escape=false)
        lag_time_name = Symbolics.tosymbol(lag_time, escape=false)

        new(
            NamedTuple{(flux_name_1,)}([fluxes[1]]),
            NamedTuple{(flux_name_2,)}([fluxes[2]]),
            NamedTuple{(lag_time_name,)}([lag_time]),
            inner_flux_func,
        )
    end
end

"""
$(TYPEDEF)
A hydrological flux calculated via a neural network (based on `Lux.jl`)
# Fields
$(FIELDS)
# Example
```
et_ann = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)
etnn_flux = NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, param_names=:etnn, chain=et_ann)
```
"""
struct NeuralFlux <: AbstractNeuralFlux
    "input variables information"
    input_info::NamedTuple
    "output variables information"
    output_info::NamedTuple
    "parameters information"
    param_info::NamedTuple
    "nn input and output information"
    nn_info::NamedTuple
    "predict function created by the chain"
    inner_func::Function
    "flux expression"
    flux_expr

    function NeuralFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        chain::Pair{Symbol,<:Lux.AbstractExplicitContainerLayer},
    )
        input_names, output_names = flux_names[1], flux_names[2]
        #* 根据输入输出参数构建 ModelingToolkit 系统
        input_vars = [first(@variables $input_name(t) = 0.0) for input_name in input_names]
        output_vars = [first(@variables $output_name(t) = 0.0) for output_name in output_names]

        return NeuralFlux(input_vars => output_vars, chain)
    end

    function NeuralFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        chain::Pair{Symbol,<:Lux.AbstractExplicitContainerLayer},
    )
        input_vars, output_vars = fluxes[1], fluxes[2]
        chain_name, chain_model = chain[1], chain[2]
        init_params = Lux.initialparameters(StableRNG(42), chain_model)

        chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(ComponentVector(init_params)))
        chain_params_type_var = first(@parameters ptype::typeof(typeof(init_params)) = typeof(init_params) [tunable = false])
        lazyconvert_params = Symbolics.array_term(convert, chain_params_type_var, chain_params, size=size(chain_params))

        input_names = Symbolics.tosymbol.(input_vars, escape=false)
        output_names = Symbolics.tosymbol.(output_vars, escape=false)

        nn_input_name = Symbol(chain_name, :_input)
        nn_output_name = Symbol(chain_name, :_output)
        nn_input = first(@variables $(nn_input_name)(t)[1:length(input_names)])
        nn_output = first(@variables $(nn_output_name)(t)[1:length(output_names)])

        #* 根据chain构建计算函数
        func = (x, p) -> LuxCore.stateless_apply(chain_model, x, p)
        flux_expr = func(nn_input, lazyconvert_params)

        function inner_flux_func(input::AbstractVector, params::AbstractVector)
            func(input, params[chain_name])
        end

        function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
            output = func(input', params[1])
            output'
        end

        new(
            NamedTuple{Tuple(input_names)}(input_vars),
            NamedTuple{Tuple(output_names)}(output_vars),
            NamedTuple{(chain_name,:ptype)}([chain_params,chain_params_type_var]),
            (input=nn_input, output=nn_output),
            inner_flux_func,
            flux_expr
        )
    end
end