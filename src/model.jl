"""
    HydroModel <: AbstractModel

Represents a hydrological model composed of multiple components.

# Fields
- `infos::NamedTuple`: Contains metadata about the model, including name, input variables, all variables, output variables, state variables, and neural network variables.
- `components::Vector{<:AbstractComponent}`: A vector of hydrological computation elements (components) that make up the model.
- `varindices::Vector`: A vector of indices for each component's input, used to map overall model inputs to component-specific inputs.

# Constructor
    HydroModel(name; components::Vector{<:AbstractComponent})

Constructs a HydroModel with the given name and components.

# Description
HydroModel is a structure that encapsulates a complete hydrological model. It manages multiple hydrological components, 
handles the flow of data between these components, and provides methods for running simulations.

The model automatically determines the connections between components based on their input and output variables. 
It also keeps track of all variables in the system, including inputs, outputs, states, and any neural network parameters.

When called as a function, the HydroModel applies its components in sequence, passing the outputs of earlier components 
as inputs to later ones, effectively simulating the hydrological system over time.

Each component's kwargs may be different, include solver, interp
"""
struct HydroModel{M<:HydroMeta,C<:AbstractComponent,VI<:AbstractVector{<:AbstractVector{<:Integer}},VN<:AbstractVector{<:Symbol}} <: AbstractModel
    "meta data of hydrological model"
    meta::M
    "hydrological computation elements"
    components::Vector{C}
    "input variables index for each components"
    varindices::VI
    "all variables names"
    varnames::VN

    function HydroModel(;
        name::Symbol,
        components::Vector{C},
        sort_components::Bool=false
    ) where {C<:AbstractComponent}
        components = sort_components ? sort_components(components) : components
        input_names, output_names, state_names = get_var_names(components)
        nn_names = reduce(union, get_nn_names.(components))
        param_names = reduce(union, get_param_names.(components))
        var_names = input_names
        input_idx = Vector{Int}[]
        for component in components
            tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, var_names), get_input_names(component))
            var_names = reduce(vcat, [var_names, get_state_names(component), get_output_names(component)])
            push!(input_idx, tmp_input_idx)
        end
        model_meta = HydroMeta(name, input_names, output_names, param_names, state_names, nn_names)
        new{typeof(model_meta),C,typeof(input_idx),typeof(var_names)}(
            model_meta,
            components,
            input_idx,
            var_names,
        )
    end
end

# 求解并计算
function (model::HydroModel)(
    input::AbstractArray{T,2}, 
    pas::ComponentVector;
    config::Union{NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
) where {T<:Number}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    fluxes = input
    for (comp_, idx, config_) in zip(model.components, model.varindices, comp_configs)
        tmp_fluxes = comp_(fluxes[idx, :], pas; config=config_, convert_to_ntp=false)
        fluxes = cat(fluxes, tmp_fluxes, dims=1)
    end
    return fluxes
end

function (model::HydroModel)(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
) where {T<:Number}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    
    fluxes = input
    for (comp_, idx_, config_) in zip(model.components, model.varindices, comp_configs)
        tmp_fluxes = comp_(fluxes[idx_, :, :], pas; config=config_, convert_to_ntp=false)
        fluxes = cat(fluxes, tmp_fluxes, dims=1)
    end
    return fluxes
end
