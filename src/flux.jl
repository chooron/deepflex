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
    input_names::Union{Symbol,Vector{Symbol}}
    "output names of the simple hydrological flux"
    output_names::Union{Symbol,Vector{Symbol}}
    "parameter names of the simple hydrological flux"
    param_names::Vector{Symbol}
    """
    functions used to implement hydrological flux computation,
    It requires the input format to be (i::NamedTuple, p::NamedTuple; kw....)
    """
    func::Function
    "smooth function used for flux func, where there exists ifelse condition"
    smooth_func::Function

    function SimpleFlux(
        input_names::Union{Symbol,Vector{Symbol}},
        output_names::Union{Symbol,Vector{Symbol}};
        param_names::Vector{Symbol},
        func::Function,
        smooth_func::Function=step_func,
    )
        return new(
            input_names,
            output_names,
            param_names,
            func,
            smooth_func
        )
    end
end

StdMeanNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:mean_, input_names), Symbol(:std_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple; kw...) ->
        @.((i[input_names] - p[param_names[1]]) / p[param_names[2]])
)


MinMaxNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:max_, input_names), Symbol(:min_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple; kw...) ->
        @.((i[input_names] - p[param_names[2]]) / (p[param_names[1]] - p[param_names[2]]))
)

TranparentFlux(
    old_names::Symbol,
    new_names::Symbol,
) = SimpleFlux(
    old_names,
    new_names,
    param_names=Symbol[],
    func=(i::NamedTuple, p::NamedTuple; kw...) -> i[old_names]
)

TranparentFlux(
    old_names::Vector{Symbol},
    new_names::Vector{Symbol},
) = SimpleFlux(
    old_names,
    new_names,
    param_names=Symbol[],
    func=(i::NamedTuple, p::NamedTuple; kw...) -> [i[nm] for nm in old_names]
)

function (flux::AbstractSimpleFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    tmp_output = flux.func(tmp_input, tmp_params, smooth_func=flux.smooth_func)
    tmp_output_names = get_output_names(flux)
    process_output(tmp_output, ifelse(length(tmp_output_names) > 1, tmp_output_names, first(tmp_output_names)))
end

"""
$(TYPEDEF)
A variable hydrological state-like flux that is determined by input and output fluxes
# Fields
$(FIELDS)
# Example
```
snowwater_flux = StateFlux([:snowfall], [:melt], :snowwater)
```
"""
mutable struct StateFlux <: AbstractStateFlux
    "input hydrological flux name for state flux"
    influx_names::Vector{Symbol}
    "output hydrological flux name for state flux"
    outflux_names::Vector{Symbol}
    "name of the state flux"
    state_names::Symbol
    "parameter names of the state flux"
    param_names::Vector{Symbol}
    """
    Function used to calculate hydrological state fluxes, fixed in most cases.
    If there is a need for modification, it is generally recommended to modify the calculation formula of the input and output flux.
    """
    func::Function

    function StateFlux(
        influx_names::Vector{Symbol},
        outflux_names::Vector{Symbol},
        state_names::Symbol;
        param_names::Vector{Symbol}=Symbol[],
        func::Function=(i::NamedTuple, p::NamedTuple; kw...) ->
            sum([i[nm] for nm in kw[:influx_names]]) - sum([i[nm] for nm in kw[:outflux_names]])
    )
        return new(
            influx_names,
            outflux_names,
            state_names,
            param_names,
            func
        )
    end
end

function (flux::AbstractStateFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    tmp_output = flux.func(tmp_input, tmp_params, influx_names=flux.influx_names, outflux_names=flux.outflux_names)
    tmp_output_names = get_output_names(flux)
    process_output(tmp_output, ifelse(length(tmp_output_names) > 1, tmp_output_names, first(tmp_output_names)))
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
    lag_func::Function
    "smooth function used for function"
    smooth_func::Function

    function LagFlux(
        input_names::Symbol,
        output_names::Symbol;
        lag_time::Symbol,
        lag_func::Function,
        smooth_func::Function=step_func,
    )
        new(
            input_names,
            output_names,
            lag_time,
            lag_func,
            smooth_func
        )
    end
end

function solve_lag_weights(flux::LagFlux; input::NamedTuple, params::Union{ComponentVector,NamedTuple})
    #* 首先将lagflux转换为discrete problem
    function discrete_prob!(du, u, p, t)
        tmp_u = input[Int(t)] .* p + u
        tmp_u = circshift(tmp_u, -1)
        tmp_u[end] = 0
        du[1:end] = tmp_u
    end
    lag_time = params[first(get_param_names(flux))]
    delta_t = 1.0
    ts = 1:(ceil(lag_time / delta_t)|>Int)
    #* 将weight作为param输入到prob中
    lag_weights = [flux.lag_func(t, lag_time, smooth_func=flux.smooth_func) for t in ts]
    lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
    prob = DiscreteProblem(discrete_prob!, zeros(eltype(eltype(input)), length(ts)), (1.0, length(input)), lag_weights)
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    sol_u = hcat((sol.u)...)
    #* 得到权重计算结果
    sol_u[1, :]
end

function (flux::LagFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    lag_weight = solve_lag_weights(flux, input=tmp_input, params=tmp_params)
    tmp_output_names = get_output_names(flux)
    process_output(lag_weight .* tmp_input[first(get_input_names(flux))],
        ifelse(length(tmp_output_names) > 1, tmp_output_names, first(tmp_output_names)))
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
    output_names::Union{Vector{Symbol},Symbol}
    "chain name of the chain inner the neural flux"
    chain_name::Symbol
    "predict function created by the chain"
    func::Function
    "prebuild neural system created by the chain(based on `ModelingToolkitNeuralNets.jl`)"
    nn::ODESystem
end

function NeuralFlux(
    input_names::Vector{Symbol},
    output_names::Union{Symbol,Vector{Symbol}};
    chain_name::Symbol,
    chain::Lux.AbstractExplicitLayer,
    seed::Int=42,
)
    func = (x, p) -> LuxCore.stateless_apply(chain, x, p)
    nn = NeuralNetworkBlock(length(input_names), length(output_names); chain=chain, rng=StableRNG(seed), name=chain_name)
    NeuralFlux(input_names, output_names, chain_name, func, nn)
end

function (flux::NeuralFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    x = hcat([input[nm] for nm in flux.input_names]...)'
    y_pred = flux.func(x, params[flux.chain_name])
    process_output(y_pred, get_output_names(flux))
end

export SimpleFlux, StateFlux, LagFlux, NeuralFlux, StdMeanNormFlux, MinMaxNormFlux, TranparentFlux