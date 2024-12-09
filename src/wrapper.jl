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
struct RecordComponentState{C<:AbstractComponent,T<:Number,M<:HydroMeta} <: AbstractHydroWrapper
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
    reset_state = get(kwargs, :reset_state, false)
    if reset_state
        wrapper.states = pas[:initstates]
    end
    pas_replace_states = update_ca(pas, ComponentVector(initstates=wrapper.states))
    output = wrapper.component(input, pas_replace_states; kwargs...)
    state_names = get_state_names(wrapper.component)
    if get(kwargs, :convert_to_ntp, false)
        component_nodes = get_node_names(wrapper.component)
        if isnothing(component_nodes)
            states_ca = ComponentVector(NamedTuple{Tuple(state_names)}([output[nm][end] for nm in state_names]))
        else
            states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
                NamedTuple{Tuple(state_names)}([output[node][nm][end] for nm in state_names])
            end))
        end
    else
        if output isa AbstractArray{<:Number,2}
            var_names = wrapper.component.varnames
            state_indices = map(state_names) do sname
                findfirst(vname -> vname == sname, var_names)
            end
            latest_states = output[state_indices, end]
            states_ca = ComponentVector(NamedTuple{Tuple(state_names)}(latest_states))
        elseif output isa AbstractArray{<:Number,3}
            latest_states = output[state_indices, :, end]
            component_nodes = get_node_names(wrapper.component)
            states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
                NamedTuple{Tuple(state_names)}(latest_states[:, node])
            end
            ))
        end
    end
    @reset wrapper.states = states_ca
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
    outlet::Integer
    meta::M

    function ComputeComponentOutlet(component::C, outlet::I) where {C<:AbstractComponent,I<:Integer}
        new{C,typeof(component.meta)}(component, outlet, component.meta)
    end
end

function (wrapper::ComputeComponentOutlet)(input::AbstractVector{<:NamedTuple}, pas::ComponentVector; kwargs...)
    convert_to_ntp = get(kwargs, :convert_to_ntp, true)
    output = wrapper.component(input, pas; kwargs...)
    return convert_to_ntp ? output[wrapper.outlet] : output[:, wrapper.outlet, :]
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

function (wrapper::WeightSumComponentOutlet)(input::Union{AbstractArray{<:Number,3},AbstractVector{<:NamedTuple}}, pas::ComponentVector; kwargs...)
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

struct NeuralWrapper{N,F,M} <: AbstractNeuralWrapper
    model::N
    func::F
    meta::M

    function NeuralWrapper(fluxes::Pair, model::N; name::Union{Nothing,Symbol}=nothing, rng::R=StableRNG(42)) where {N,R}
        #* assert the chain has a name
        # @assert model.name isa Symbol "Neural network chain should have a name with Symbol type"
        #* extract the parameter and state
        ps, st = Lux.setup(rng, model)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func(x, p) = LuxCore.apply(model, x, ComponentVector(p, ps_axes), st)[1]
        wrapper_name = name === nothing ? Symbol(model.name, :_wrapper) : name
        meta = HydroMeta(name=wrapper_name, inputs=fluxes.first, outputs=fluxes.second) # , nn_names=[model.name]
        new{N,typeof(nn_func),typeof(meta)}(model, nn_func, meta)
    end
end

function (wrapper::NeuralWrapper)(input::Array{T,2}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nn][wrapper.meta.name]
    output = wrapper.func(input, nn_ps)
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)
    return convert_to_ntp ? NamedTuple{Tuple(get_output_names(wrapper))}(eachslice(output, dims=1)) : output
end

function (wrapper::NeuralWrapper)(input::Array{T,3}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nn][wrapper.meta.name]
    output = wrapper.func.(Matrix.(eachslice(input, dims=2)), Ref(nn_ps))
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)
    if convert_to_ntp
        return [NamedTuple{Tuple(get_output_names(wrapper))}(eachslice(output_, dims=1)) for output_ in eachslice(output, dims=2)]
    else
        return permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (1,3,2))
    end
end