function (flux::AbstractFlux)(input::Vector, params::Vector)
    flux.inner_func(input, params)
end

function (flux::AbstractFlux)(input::Matrix, params::Vector)
    flux.inner_func(input, params)
end
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
        input_names::Union{Symbol,Vector{Symbol}},
        output_names::Union{Symbol,Vector{Symbol}};
        param_names::Vector{Symbol}=Symbol[],
        smooth_func::Function=step_func,
    )
        if input_names isa Symbol
            input_names = [input_names]
        end

        if output_names isa Symbol
            output_names = [output_names]
        end

        #* 根据输入输出参数名称获取对应的计算公式
        hydro_equation = HydroEquation(input_names, output_names, param_names)
        flux_exprs = expr(hydro_equation, smooth_func=smooth_func)

        #* 得到计算函数
        flux_funcs = [build_function(
            hydro_expr, hydro_equation.inputs, hydro_equation.params, expression=Val{false}
        ) for hydro_expr in flux_exprs]

        #* 得到计算公式
        flux_eqs = [output ~ hydro_expr for (output, hydro_expr) in zip(hydro_equation.outputs, flux_exprs)]

        function inner_flux_func(input::Vector, params::Vector)
            [func(input, params) for func in flux_funcs]
        end

        function inner_flux_func(input::Matrix, params::Vector)
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

    function SimpleFlux(
        inputs::Vector{Num},
        outputs::Vector{Num};
        params::Vector{Num}=Num[],
        exprs::Vector{Num}
    )
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params)

        #* 得到计算函数
        flux_funcs = [build_function(flux_expr, inputs, params, expression=Val{false}) for flux_expr in exprs]

        #* 得到计算公式
        flux_eqs = [output ~ flux_expr for (output, flux_expr) in zip(outputs, exprs)]

        function inner_flux_func(input::Vector, params::Vector)
            #* 单值计算
            [func(input, params) for func in flux_funcs]
        end

        return new(
            input_names,
            output_names,
            param_names,
            inner_flux_func,
            flux_eqs,
            exprs,
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
    influx_names::Vector{Symbol}
    "output hydrological flux name for state flux"
    outflux_names::Vector{Symbol}
    "name of the state flux"
    state_name::Symbol
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
        influx_names::Vector{Symbol},
        outflux_names::Vector{Symbol},
        state_name::Symbol;
        fluxes::Vector{<:AbstractFlux}
    )
        #* 获得elements的输入变量名
        #* 1. 这个是fluxes计算需要输入的通量
        fluxes_input_names, fluxes_output_names = get_input_output_names(fluxes)
        #* 2. 这个是state不通过fluxes计算的输入通量
        state_input_names = setdiff(vcat(influx_names, outflux_names), fluxes_output_names)
        #* 3. 合并起来就是计算state所需的最基础通量(不包括状态变量)
        union_input_names = union(fluxes_input_names, state_input_names)

        #* fluxes计算过程中需要构建的通量
        fluxes_var_names = vcat(union_input_names, fluxes_output_names)
        fluxes_param_names = get_param_names(fluxes)

        #* 构建变量
        fluxes_vars = [first(@variables $nm(t) = 0.0) for nm in fluxes_var_names]
        fluxes_params = [first(@parameters $param_name = 0.0 [tunable = true]) for param_name in fluxes_param_names]
        #* 构建state变量
        flux_state = first(@variables $state_name(t) = 0.0)
        fluxes_vars_ntp = NamedTuple{Tuple(fluxes_var_names)}(fluxes_vars)

        #* 构建函数和公式
        state_expr = sum(fluxes_vars_ntp[influx_names]) - sum(fluxes_vars_ntp[outflux_names])
        state_eq = D(flux_state) ~ state_expr

        #* 构建state计算的函数并将所有中间状态替换
        substitute_vars_dict = Dict()
        for var_nm in keys(fluxes_vars_ntp)
            for flux in fluxes
                for j in eachindex(flux.output_names)
                    if var_nm == flux.output_names[j]
                        tmp_flux_exprs = Symbolics.rhss(flux.flux_eqs)
                        substitute_vars_dict[fluxes_vars_ntp[var_nm]] = tmp_flux_exprs[j]
                    end
                end
            end
        end
        state_expr_sub = substitute(state_expr, substitute_vars_dict)
        state_func = build_function(state_expr_sub, collect(fluxes_vars_ntp[union_input_names]), fluxes_params, expression=Val{false})

        return new(
            influx_names,
            outflux_names,
            state_name,
            state_func,
            state_eq
        )
    end
end

function get_state_input_var_names(state_flux::StateFlux)

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
        input_names::Symbol,
        output_names::Symbol;
        lag_time::Symbol,
        lag_func::Function,
        smooth_func::Function=step_func,
    )
        function inner_flux_func(i::StructArray, p::NamedTuple)
            tmp_input = i[input_names]
            lag_weight = solve_lag_weights(tmp_input, p[lag_time], lag_func, smooth_func)
            i[output_names] .= lag_weight .* tmp_input
        end

        new(
            input_names,
            output_names,
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
end

function NeuralFlux(
    input_names::Union{Symbol,Vector{Symbol}},
    output_names::Union{Symbol,Vector{Symbol}};
    chain_name::Symbol,
    chain::Lux.AbstractExplicitLayer,
    seed::Int=42,
)
    if input_names isa Symbol
        input_names = [input_names]
    end
    if output_names isa Symbol
        output_names = [output_names]
    end

    #* 根据输入输出参数构建ModelingToolkit系统
    input_vars = [first(@variables $input_name(t) = 0.0) for input_name in input_names]
    output_vars = [first(@variables $output_name(t) = 0.0) for output_name in output_names]

    nn_in = RealInputArray(nin=length(input_names), name=Symbol(chain_name, :_in_sys))
    nn_out = RealOutputArray(nout=length(output_names), name=Symbol(chain_name, :_out_sys))
    nn = NeuralNetworkBlock(
        length(input_names), length(output_names);
        chain=chain, rng=StableRNG(seed), name=Symbol(chain_name, :_nn_sys)
    )

    eqs = Equation[connect(nn_in, nn.input), connect(nn_out, nn.output)]
    for i in eachindex(input_vars)
        push!(eqs, input_vars[i] ~ nn_in.u[i])
    end
    for i in eachindex(output_vars)
        push!(eqs, output_vars[i] ~ nn_out.u[i])
    end

    #* 根据chain构建计算函数
    func = (x, p) -> LuxCore.stateless_apply(chain, x, p)

    function inner_flux_func(input::NamedTuple, params::NamedTuple)
        """
        单值计算
        """
        x = hcat([input[nm] for nm in input_names]...)'
        output = func(x, params[chain_name])
        [first(output[idx, :]) for idx in 1:length(output_names)]
    end

    function inner_flux_func(input::StructArray, params::NamedTuple)
        """
        数组计算
        """
        x = hcat([getproperty(input, nm) for nm in input_names]...)'
        output = func(x, params[chain_name])
        [output[idx, :] for idx in 1:length(output_names)]
    end

    NeuralFlux(input_names, output_names, chain_name, inner_flux_func, eqs, [nn, nn_in, nn_out])
end

function extract_neuralflux_ntp(funcs::Vector{<:AbstractFlux})
    nn_flux_list = filter(flux -> flux isa AbstractNeuralFlux, funcs)
    nfunc_ntp = NamedTuple()
    flxu_chain_names = [flux.chain_name for flux in nn_flux_list]
    if length(nn_flux_list) > 0
        nfunc_ntp = NamedTuple{Tuple(flxu_chain_names)}(
            [get_input_names(flux) for flux in nn_flux_list]
        )
    end
    nfunc_ntp
end