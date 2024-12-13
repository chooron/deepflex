"""
    EstimateComponentParams <: AbstractHydroWrapper

A wrapper component that estimates parameters for a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `estfuncs::AbstractVector{<:AbstractHydroFlux}`: Vector of estimation functions
- `ptypes::AbstractVector{Symbol}`: Parameter types to be estimated
- `meta::HydroMeta`: Metadata about the component

# Constructors
    EstimateComponentParams(component::AbstractComponent, estfuncs::AbstractVector{<:AbstractHydroFlux}, ptypes::AbstractVector{Symbol}=Symbol[])

# Description
`EstimateComponentParams` wraps a hydrological component and dynamically estimates its parameters
using provided estimation functions during simulation.

The wrapper requires at least one estimation function. Each estimation function should implement:
- `get_input_names`: Returns required input variable names
- `get_param_names`: Returns parameter names it can estimate
- `get_var_names`: Returns both input and output variable names

When called, it:
1. Applies each estimation function to calculate parameters
2. Updates the parameter set with estimated values
3. Runs the wrapped component with updated parameters
4. Returns the component output

The estimation process uses the parameter types specified in `ptypes` to organize
and update the correct parameter groups in the component.
"""
struct EstimateComponentParams{C<:AbstractComponent,E<:AbstractHydroFlux,P<:Union{AbstractVector{Symbol},Nothing},M<:HydroMeta} <: AbstractHydroWrapper
    component::C
    estfuncs::AbstractVector{E}
    pkeys::P
    meta::M

    function EstimateComponentParams(
        component::C, estfuncs::AbstractVector{E}, pkeys::P=nothing
    ) where {C<:AbstractComponent,E<:AbstractHydroFlux,P<:Union{AbstractVector{Symbol},Nothing}}
        @assert length(estfuncs) != 0 "At least one estimation function is required"
        comp_meta = component.meta
        est_input_names, est_output_names = get_var_names(estfuncs)
        est_param_names = reduce(union, get_param_names.(estfuncs))
        # todo: due to ComponentVector merging limitations, arbitrary values still need to be provided for estimated parameters
        # new_param_names = setdiff(union(comp_meta.params, est_input_names, est_param_names), est_output_names)
        new_param_names = union(comp_meta.params, est_input_names, est_param_names)
        new_meta = HydroMeta(comp_meta.name, comp_meta.inputs, comp_meta.outputs, new_param_names, comp_meta.states, comp_meta.nns)
        return new{C,E,P,typeof(new_meta)}(component, estfuncs, pkeys, new_meta)
    end
end

function (wrapper::EstimateComponentParams{C,E,P,M})(input::Any, pas::ComponentVector; kwargs...) where {C,E,P<:Nothing,M}
    est_params_ntp = reduce(merge, map(wrapper.estfuncs) do estfunc
        tmp_input = collect(pas[:params][get_input_names(estfunc)])
        tmp_pas = collect(pas[:params][get_param_names(estfunc)])
        NamedTuple{Tuple(get_output_names(estfunc))}(estfunc(tmp_input, pas))
    end)
    new_pas = update_ca(pas, ComponentVector(params=est_params_ntp))
    output = wrapper.component(input, new_pas; kwargs...)
    return output
end

function (wrapper::EstimateComponentParams{C,E,P,M})(input::Any, pas::ComponentVector; kwargs...) where {C,E,P<:AbstractVector{Symbol},M}
    est_config = (ptypes=wrapper.pkeys,)
    est_params_mat = reduce(vcat, map(wrapper.estfuncs) do estfunc
        tmp_input_mat = reduce(vcat, map(wrapper.pkeys) do pkey
            tmp_input_vec = collect([pas[:params][pkey][nm] for nm in get_input_names(estfunc)])
            reshape(tmp_input_vec, 1, length(tmp_input_vec))
        end)
        tmp_input_arr = permutedims(reshape(tmp_input_mat, 1, size(tmp_input_mat)...), (3, 2, 1))
        tmp_output_mat = estfunc(tmp_input_arr, pas, config=est_config)[:, :, 1]
        tmp_output_mat
    end)
    est_param_names = reduce(union, get_output_names.(wrapper.estfuncs))
    est_params_ntp = NamedTuple{Tuple(wrapper.pkeys)}([NamedTuple{Tuple(est_param_names)}(est_params_mat[:, i]) for i in eachindex(wrapper.pkeys)])
    new_pas = update_ca(pas, ComponentVector(params=est_params_ntp))
    output = wrapper.component(input, new_pas; kwargs...)
    return output
end