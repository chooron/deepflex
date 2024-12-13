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