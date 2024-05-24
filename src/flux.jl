function (flux::AbstractFlux)(input::NamedTuple, params::NamedTuple)
    flux.inner_func(input, params)
end

function (flux::AbstractFlux)(input::StructArray, params::NamedTuple)
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
    It requires the input format to be (i::NamedTuple, p::NamedTuple; kw....)
    """
    inner_func::Function

    function SimpleFlux(
        input_names::Union{Symbol,Vector{Symbol}},
        output_names::Union{Symbol,Vector{Symbol}};
        param_names::Vector{Symbol},
        func::Function,
        smooth_func::Function=step_func,
    )
        if input_names isa Symbol
            input_names = [input_names]
        end

        if output_names isa Symbol
            output_names = [output_names]
        end

        function inner_flux_func(input::NamedTuple, params::NamedTuple)
            """
            单值计算
            """
            tmp_input = input[input_names]
            tmp_output = func(tmp_input, params[param_names], smooth_func=smooth_func)
            tmp_output
        end

        function inner_flux_func(input::StructArray, params::NamedTuple)
            """
            数组计算
            """
            tmp_input = StructArrays.components(input)[input_names]
            tmp_output = func(tmp_input, params[param_names], smooth_func=smooth_func)
            tmp_output
        end

        return new(
            input_names,
            output_names,
            param_names,
            inner_flux_func
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
mutable struct StateFlux{mtk} <: AbstractStateFlux
    "input hydrological flux name for state flux"
    influx_names::Vector{Symbol}
    "output hydrological flux name for state flux"
    outflux_names::Vector{Symbol}
    "name of the state flux"
    state_names::Symbol
    """
    Function used to calculate hydrological state fluxes, fixed in most cases.
    If there is a need for modification, it is generally recommended to modify the calculation formula of the input and output flux.
    """
    inner_func::Function

    function StateFlux(
        influx_names::Vector{Symbol},
        outflux_names::Vector{Symbol},
        state_names::Symbol;
        mtk::Bool=true
    )
        function inner_flux_func_mtk(i::NamedTuple, p::NamedTuple)
            sum([i[nm] for nm in influx_names]) - sum([i[nm] for nm in outflux_names])
        end

        return new{mtk}(
            influx_names,
            outflux_names,
            state_names,
            inner_flux_func_mtk
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
    nn::ODESystem
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
    nn = NeuralNetworkBlock(length(input_names), length(output_names); chain=chain, rng=StableRNG(seed), name=chain_name)
    NeuralFlux(input_names, output_names, chain_name, inner_flux_func, nn)
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