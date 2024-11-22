"""
    RecordComponentState <: AbstractHydroWrapper

A wrapper component that records the state values of a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `states::ComponentVector`: Current state values of the component
- `meta::HydroMeta`: Metadata about the component

# Constructors
    RecordComponentState(component::AbstractComponent, initstates::Union{Nothing,ComponentVector}=nothing)

# Description
`RecordComponentState` wraps a hydrological component and tracks its state values during simulation.
It updates and stores the latest state values after each simulation step.

The wrapper requires that the wrapped component has state variables. If no initial states are provided,
it will use default states from the component.

When called, it:
1. Updates the input parameters with current state values
2. Runs the wrapped component simulation
3. Extracts and stores the final state values
4. Returns the original component output

The state values can be accessed through the `states` field at any time.
"""
mutable struct RecordComponentState{C<:AbstractComponent,T<:Number,M<:HydroMeta} <: AbstractHydroWrapper
    component::C
    states::ComponentVector{T}
    meta::M

    function RecordComponentState(component::C, initstates::Union{Nothing,ComponentVector}=nothing) where {C<:AbstractComponent}
        @assert length(get_state_names(component)) != 0 "Component $component must have state flux"
        initstates = isnothing(initstates) ? get_default_states(component) : initstates
        new{C,eltype(initstates),typeof(component.meta)}(component, initstates, component.meta)
    end
end

function (wrapper::RecordComponentState{C,T,M})(input::Any, pas::ComponentVector; kwargs...) where {C,T,M}
    pas_replace_states = update_ca(pas, wrapper.states)
    output = wrapper.component(input, pas_replace_states; kwargs...)
    state_names = get_state_names(wrapper.component)
    if get(kwargs, :convert_to_ntp, false)
        component_nodes = get_node_names(wapper.component)
        if isnothing(component_nodes)
            states_ca = ComponentVector(NamedTuple{Tuple(state_names)}([output[nm][end] for nm in state_names]))
        else
            states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
                NamedTuple{Tuple(state_names)}([output[node][nm][end] for nm in state_names])
            end))
        end
    else
        if output isa AbstractArray{<:Number,2}
            var_names = wapper.component.varnames
            state_indices = map(state_names) do sname
                findfirst(vname -> vname == sname, var_names)
            end
            latest_states = output[state_indices, end]
            states_ca = ComponentVector(NamedTuple{Tuple(state_names)}(latest_states))
        elseif output isa AbstractArray{<:Number,3}
            latest_states = output[state_indices, :, end]
            component_nodes = get_node_names(wapper.component)
            states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
                NamedTuple{Tuple(state_names)}(latest_states[:, node])
            end
            ))
        end
    end
    wapper.states = states_ca
    return output
end

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
struct EstimateComponentParams{C<:AbstractComponent,E<:AbstractHydroFlux,M<:HydroMeta} <: AbstractHydroWrapper
    component::C
    estfuncs::AbstractVector{E}
    ptypes::AbstractVector{Symbol}
    meta::M

    function EstimateComponentParams(component::C, estfuncs::AbstractVector{E}, ptypes::AbstractVector{Symbol}=Symbol[]) where {C<:AbstractComponent,E<:AbstractHydroFlux}
        @assert length(estfuncs) != 0 "At least one estimation function is required"
        comp_meta = component.meta
        est_input_names, est_output_names = get_var_names(estfuncs)
        est_param_names = reduce(union, get_param_names.(estfuncs))
        new_param_names = setdiff(union(comp_meta.param_names, est_input_names, est_param_names), est_output_names)
        new_meta = HydroMeta(comp_meta.name, comp_meta.input_names, comp_meta.output_names, new_param_names, comp_meta.state_names, comp_meta.nn_names)
        return new{C,E,typeof(new_meta)}(component, estfuncs, ptypes, new_meta)
    end
end

function (wrapper::EstimateComponentParams)(input::Any, pas::ComponentVector; kwargs...)
    est_params_arr = reduce(hcat, map(wrapper.estfuncs) do estfunc
        tmp_input = [pas[:params][ptype][get_input_names(estfunc)] for ptype in wrapper.ptypes]
        tmp_pas = [pas[:params][ptype][get_param_names(estfunc)] for ptype in wrapper.ptypes]
        estfunc.(tmp_input, tmp_pas; convert_to_ntp=false)
    end)
    total_est_param_names = reduce(union, get_param_names.(wrapper.estfuncs))
    est_params_ca = NamedTuple{Tuple(wrapper.ptypes)}([NamedTuple{Tuple(total_est_param_names)}(est_params_arr[:, i]) for i in eachindex(wrapper.ptypes)])
    new_pas = update_ca(pas, est_params_ca)
    output = wrapper.component(input, new_pas; kwargs...)
    return output
end

"""
    ComputeComponentOutlet <: AbstractHydroWrapper

A wrapper component that computes the outlet values of a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `outlet::Symbol`: The name of the outlet node to compute values for
- `meta::HydroMeta`: Metadata about the component

# Constructors
    ComputeComponentOutlet(component::AbstractComponent, outlet::Symbol, meta::HydroMeta)

# Description
`ComputeComponentOutlet` wraps a hydrological component and extracts the values at a specific outlet node
during simulation.

The wrapper requires that the outlet name exists in the component's node names. When called, it:
1. Runs the wrapped component simulation
2. Extracts the values at the specified outlet node
3. Returns just the outlet values, either as a named tuple or array depending on configuration

This is useful for focusing analysis on a particular outlet point in a hydrological network.
"""
struct ComputeComponentOutlet{C<:AbstractComponent,M<:HydroMeta} <: AbstractHydroWrapper
    component::C
    outlet::Symbol
    meta::M

    function ComputeComponentOutlet(component::C, outlet::S) where {C<:AbstractComponent,S<:Symbol}
        @assert outlet in component.hrunames "Outlet $outlet is not a valid output name for component $component"
        new{C,typeof(meta)}(component, outlet, meta)
    end
end

function (wrapper::ComputeComponentOutlet)(input::AbstractVector{<:NamedTuple}, pas::ComponentVector; kwargs...)
    convert_to_ntp = get(kwargs, :convert_to_ntp, true)
    @assert wrapper.outlet in wrapper.component.varnames "Outlet $wrapper.outlet is not a valid output name for component $wrapper.component"
    outlet_idx = findfirst(nm -> nm == wrapper.outlet, wrapper.component.varnames)
    output = wrapper.component(input, pas; kwargs...)
    return convert_to_ntp ? output[outlet_idx] : output[:, outlet_idx, :]
end

"""
    ComputeComponentOutlet <: AbstractHydroWrapper

A wrapper component that computes the outlet values of a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `outlet::Symbol`: The name of the outlet node to compute values for
- `meta::HydroMeta`: Metadata about the component

# Constructor
    ComputeComponentOutlet(component::AbstractComponent, outlet::Symbol)

# Description
`ComputeComponentOutlet` wraps a hydrological component and extracts the values at a specific outlet node
during simulation.

The wrapper requires that the outlet name exists in the component's node names. When called, it:
1. Runs the wrapped component simulation
2. Extracts the values at the specified outlet node
3. Returns just the outlet values, either as a named tuple or array depending on configuration

This is useful for focusing analysis on a particular outlet point in a hydrological network.
"""
struct WeightSumComponentOutlet{C<:AbstractComponent,M<:HydroMeta,T<:Number} <: AbstractHydroWrapper
    component::C
    vars::AbstractVector{T}
    meta::M
    weight_names::Symbol

    function WeightSumComponentOutlet(component::C, vars::AbstractVector{T}; weight_names::Symbol=:weight) where {C<:AbstractComponent,T<:Num}
        new_params = vcat(component.meta.param_names, [weight_names])
        var_names = Symbolics.tosymbol.(vars)
        new_meta = HydroMeta(component.meta.name, component.meta.input_names, new_params, var_names, component.meta.state_names, component.meta.nn_names)
        new{C,typeof(meta),T}(component, vars, new_meta, weight_names)
    end
end

function (wrapper::WeightSumComponentOutlet)(input::Union{AbstractArray{T,3},AbstractVector{<:NamedTuple}}, pas::ComponentVector; kwargs...) where {T}
    convert_to_ntp = get(kwargs, :convert_to_ntp, true)
    output = wrapper.component(input, pas; kwargs...)
    var_names = Symbolics.tosymbol.(wrapper.vars)
    weight_vec = [pas[:params][ptype][wrapper.weight_names] for ptype in wrapper.ptypes]
    if convert_to_ntp
        var_vec = map(var_names) do var_nm
            sum(reduce(hcat, [out[var_nm] .* weight_vec[i] for (i, out) in enumerate(weight_vec)]))
        end
        return NamedTuple{Tuple(var_names)}(var_vec)
    else
        var_idx = [findfirst(nm -> nm == var_nm, wrapper.component.varnames) for var_nm in var_names]
        var_vec = map(var_idx) do idx
            sum([output[idx, i, :] .* weight_vec[i] for i in eachindex(weight_vec)])
        end
        return reduce(hcat, var_vec)
    end
end