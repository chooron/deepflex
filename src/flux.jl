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
    "input names of the simple hydrological flux"
    input_names::Vector{Symbol}
    "output names of the simple hydrological flux"
    output_names::Vector{Symbol}
    "parameter names of the simple hydrological flux"
    param_names::Vector{Symbol}
    """
    functions used to implement hydrological flux computation,
    It requires the input format to be (i::Vector, p::Vector)
    """
    inner_func::Function
    """
    flux equations
    """
    flux_eqs::Vector{Equation}

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
            params = [first(@parameters $var = 0.0) for var in output_names]
            flux_eqs = [output ~ flux_func(inputs, params) for (output, flux_func) in zip(outputs, flux_funcs)]
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
            flux_eqs = [output ~ hydro_expr for (output, hydro_expr) in zip(hydro_equation.outputs, flux_flux_exprs)]
        end

        function inner_flux_func(input::AbstractArray, params::AbstractArray)
            [func(input, params) for func in flux_funcs]
        end

        # function inner_flux_func(input::AbstractMatrix, params::AbstractVector, timeidx::Vector)
        #     #! 这个地方有问题
        #     [[func(input[i, :], params) for i in timeidx] for func in flux_funcs]
        # end

        # function inner_flux_func(input::NamedTuple, params::NamedTuple)
        #     [[func([input[nm] for nm in input_names], collect(params[param_names]))] for func in flux_funcs]
        # end

        return new(
            input_names,
            output_names,
            param_names,
            inner_flux_func,
            flux_eqs,
        )
    end

    function SimpleFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        flux_exprs::Vector{Num}
    )
        inputs, outputs = fluxes[1], fluxes[2]

        if inputs isa Num
            inputs = [inputs]
        end

        if outputs isa Num
            outputs = [outputs]
        end

        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params)

        #* 得到计算函数
        flux_funcs = [build_function(flux_expr, inputs, params, expression=Val{false}) for flux_expr in flux_exprs]

        #* 得到计算公式
        flux_eqs = [output ~ flux_expr for (output, flux_expr) in zip(outputs, flux_exprs)]

        function inner_flux_func(input::AbstractArray, params::AbstractArray)
            [func(input, params) for func in flux_funcs]
        end

        function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
            [[func(input[i, :], params) for i in 1:size(input)[1]] for func in flux_funcs]
        end

        return new(
            input_names,
            output_names,
            param_names,
            inner_flux_func,
            flux_eqs,
        )
    end
end

"""
$(TYPEDEF)
A variable hydrological state-like flux that is determined by input and output fluxes
# Fields
$(FIELDS)
# Example
```
snowwater_flux = StateFlux([:snowfall], [:melt], :snowwater, [])
```
"""
struct StateFlux <: AbstractStateFlux
    "input hydrological flux name for state flux"
    input_names::Vector{Symbol}
    "output hydrological flux name for state flux"
    output_names::Vector{Symbol}
    """
    Function used to calculate hydrological state fluxes, fixed in most cases.
    If there is a need for modification, it is generally recommended to modify the calculation formula of the input and output flux.
    """
    inner_func::Function
    """
    Equation used to calculate hydrological state fluxes, fixed in most cases.
    """
    state_eq::Equation

    function StateFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        state::Num;
        funcs::Vector{<:AbstractFlux}
    )
        influxes, outfluxes = fluxes[1], fluxes[2]
        #* 转换为Symbol
        influx_names = Symbolics.tosymbol.(influxes, escape=false)
        outflux_names = Symbolics.tosymbol.(outfluxes, escape=false)
        state_input_names = Symbol.(union(vcat(influx_names, outflux_names)))
        state_name = Symbolics.tosymbol(state, escape=false)

        #* 构建函数和公式
        state_expr = sum(influxes) - sum(outfluxes)
        state_eq = D(state) ~ state_expr

        #* 获取state flux name信息
        funcs_input_names, funcs_output_names = get_input_output_names(funcs)
        union_input_names = union(funcs_input_names, setdiff(state_input_names, funcs_output_names))
        funcs_var_names = vcat(union_input_names, funcs_output_names)
        funcs_param_names = get_param_names(funcs)

        #* 构建计算函数
        state_func = build_state_func(funcs, state_expr, union_input_names, funcs_var_names, funcs_param_names)

        return new(
            union_input_names,
            [state_name],
            state_func,
            state_eq
        )
    end

    function StateFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        state_name::Symbol;
        funcs::Vector{<:AbstractFlux}
    )
        influx_names, outflux_names = flux_names[1], flux_names[2]
        state_input_names = vcat(influx_names, outflux_names)
        #* 构建函数和公式
        influxes = [first(@variables $nm(t) = 0.0) for nm in influx_names]
        outfluxes = [first(@variables $nm(t) = 0.0) for nm in outflux_names]
        state = first(@variables $state_name(t) = 0.0)
        state_expr = sum(influxes) - sum(outfluxes)
        state_eq = D(state) ~ state_expr

        #* 获取state flux name信息
        funcs_input_names, funcs_output_names = get_input_output_names(funcs)
        union_input_names = union(funcs_input_names, setdiff(state_input_names, funcs_output_names))
        funcs_var_names = vcat(union_input_names, funcs_output_names)
        funcs_param_names = get_param_names(funcs)

        #* 构建计算函数
        state_func = build_state_func(funcs, state_expr, union_input_names, funcs_var_names, funcs_param_names)

        return new(
            union_input_names,
            [state_name],
            state_func,
            state_eq
        )
    end

    function StateFlux(
        state_names::Pair{Symbol,Symbol};
        funcs::Vector{<:AbstractFlux}
    )
        new_state_name, ori_state_name = state_names[1], state_names[2]
        #* 构建函数和公式
        ori_state = first(@variables $ori_state_name(t) = 0.0)
        new_state = first(@variables $new_state_name(t) = 0.0)
        state_expr = new_state - ori_state
        state_eq = D(ori_state) ~ state_expr

        #* 构建计算函数
        fluxes_input_names = get_input_names(fluxes)
        new_state_func = build_new_state_func(funcs, state_expr)

        return new(
            setdiff(fluxes_input_names, [ori_state_name]),
            [ori_state_name],
            new_state_func,
            state_eq
        )
    end

    function StateFlux(
        states::Pair{Num,Num};
        funcs::Vector{<:AbstractFlux}
    )
        #* 构建函数和公式
        new_state, ori_state = states[1], states[2]
        state_expr = new_state - ori_state
        state_eq = D(ori_state) ~ state_expr

        ori_state_name = Symbolics.tosymbol(ori_state, escape=false)
        new_state_name = Symbolics.tosymbol(new_state, escape=false)

        #* 构建计算函数
        fluxes_input_names = get_input_names(funcs)
        new_state_func = build_new_state_func(funcs, state_expr)

        return new(
            fluxes_input_names,
            [ori_state_name],
            new_state_func,
            state_eq
        )
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
    "input name of the hydrograph unit flux, only supports a single name"
    input_names::Symbol
    "output name of the hydrograph unit flux, only supports a single name"
    output_names::Symbol
    "parameter name used to represent the unit line parameter (lag time)"
    lag_time::Symbol
    "function used to represent hydrograph unit"
    inner_func::Function

    function LagFlux(
        flux_name::Symbol,
        lag_time::Symbol,
        lag_func::Function;
        kwargs...,
    )
        function inner_flux_func(i::AbstractMatrix, p::AbstractVector)
            lag_flux = solve_lag_flux(i[:, 1], p[1], lag_func, kwargs...)
            [lag_flux]
        end

        new(
            flux_name,
            Symbol(flux_name, :_lag),
            lag_time,
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
    "input names of the neural flux"
    input_names::Vector{Symbol}
    "output names of the neural flux"
    output_names::Vector{Symbol}
    "chain name of the chain inner the neural flux"
    chain_name::Symbol
    "predict function created by the chain"
    inner_func::Function
    "prebuild neural system created by the chain(based on `ModelingToolkitNeuralNets.jl`)"
    flux_eqs::Vector{Equation}
    "neural network sub systems"
    sub_sys::Vector{ODESystem}

    function NeuralFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        chain::Pair{Symbol,Lux.AbstractExplicitLayer},
        seed::Int=42,
    )
        input_names, output_names = flux_names[1], flux_names[2]
        if input_names isa Symbol
            input_names = [input_names]
        end
        if output_names isa Symbol
            output_names = [output_names]
        end

        #* 根据输入输出参数构建ModelingToolkit系统
        input_vars = [first(@variables $input_name(t) = 0.0) for input_name in input_names]
        output_vars = [first(@variables $output_name(t) = 0.0) for output_name in output_names]

        return NeuralFlux(input_vars => output_vars, chain_name, chain, seed)
    end

    function NeuralFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        chain::Pair{Symbol,Lux.AbstractExplicitLayer},
        seed::Int=42,
    )
        input_vars, output_vars = fluxes[1], fluxes[2]
        chain_name, chain_model = chain[1], chain[2]

        if input_vars isa Num
            input_vars = [input_vars]
        end
        if output_vars isa Num
            output_vars = [output_vars]
        end

        input_names = Symbolics.tosymbol.(input_vars, escape=false)
        output_names = Symbolics.tosymbol.(output_vars, escape=false)

        nn_in = RealInputArray(nin=length(input_names), name=Symbol(chain_name, :_in_sys))
        nn_out = RealOutputArray(nout=length(output_names), name=Symbol(chain_name, :_out_sys))
        nn = NeuralNetworkBlock(
            length(input_names), length(output_names);
            chain=chain_model, rng=StableRNG(seed), name=Symbol(chain_name, :_nn_sys)
        )

        eqs = Equation[connect(nn_in, nn.input), connect(nn_out, nn.output)]
        for i in eachindex(input_vars)
            push!(eqs, input_vars[i] ~ nn_in.u[i])
        end
        for i in eachindex(output_vars)
            push!(eqs, output_vars[i] ~ nn_out.u[i])
        end

        #* 根据chain构建计算函数
        func = (x, p) -> LuxCore.stateless_apply(chain_model, x, p)

        function inner_flux_func(input::AbstractVector, params::AbstractVector)
            x = hcat(input...)'
            output = func(x, params[chain_name])
            [first(output[idx, :]) for idx in 1:length(output_names)]
        end

        new(
            input_names,
            output_names,
            chain_name,
            inner_flux_func,
            eqs,
            [nn, nn_in, nn_out]
        )
    end
end