"""
    HydroMeta

A structure that stores metadata about hydrological components, including variable names and component information.

# Fields
- `inputs::Vector{Symbol}`: Names of input variables
- `outputs::Vector{Symbol}`: Names of output variables
- `states::Vector{Symbol}`: Names of state variables
- `params::Vector{Symbol}`: Names of parameters
- `nns::Vector{Symbol}`: Names of neural network components
- `name::Symbol`: Component name
"""
struct HydroMeta
    name::Symbol
    inputs::NamedTuple
    outputs::NamedTuple
    params::NamedTuple
    states::NamedTuple
    nns::NamedTuple

    function HydroMeta(
        name::Symbol, input_names::Vector{Symbol}=Symbol[], output_names::Vector{Symbol}=Symbol[],
        param_names::Vector{Symbol}=Symbol[], state_names::Vector{Symbol}=Symbol[], nn_names::Vector{Symbol}=Symbol[],
    )
        return new(name, input_names, output_names, param_names, state_names, nn_names)
    end

    function HydroMeta(;
        name::Symbol, inputs::Vector{Num}=Num[], outputs::Vector{Num}=Num[],
        params::Vector{Num}=Num[], states::Vector{Num}=Num[], nn_names::Vector{Symbol}=Symbol[],
    )
        input_names = isempty(inputs) ? Symbol[] : tosymbol.(inputs, escape=false)
        output_names = isempty(outputs) ? Symbol[] : tosymbol.(outputs, escape=false)
        param_names = isempty(params) ? Symbol[] : tosymbol.(params, escape=false)
        state_names = isempty(states) ? Symbol[] : tosymbol.(states, escape=false)
        return new(name, input_names, output_names, param_names, state_names, nn_names)
    end
end


"""
    get_names(comp::AbstractComponent)

Get the name of a component.
Returns a symbol representing the component name.
"""
get_name(comp::AbstractComponent) = comp.meta.name

"""
    get_input_names(comp::AbstractComponent)
    get_input_names(comps::Vector{<:AbstractComponent})

Get the names of input variables from a component or vector of components.
Returns a vector of symbols representing input variable names.
"""
get_input_names(comp::AbstractComponent) = comp.meta.inputs
get_input_names(comps::Vector{<:AbstractComponent}) = reduce(union, get_input_names.(comps))

"""
    get_output_names(comp::AbstractComponent)
    get_output_names(comps::Vector{<:AbstractComponent})

Get the names of output variables from a component or vector of components.
Returns a vector of symbols representing output variable names.
"""
get_output_names(comp::AbstractComponent) = comp.meta.outputs
get_output_names(comps::Vector{<:AbstractComponent}) = reduce(union, get_output_names.(comps))

"""
    get_state_names(comp::AbstractComponent)
    get_state_names(comps::Vector{<:AbstractComponent})

Get the names of state variables from a component or vector of components.
Returns a vector of symbols representing state variable names.
"""
get_state_names(comp::AbstractComponent) = comp.meta.states
get_state_names(comps::Vector{<:AbstractComponent}) = reduce(union, get_state_names.(comps))

"""
    get_param_names(comp::AbstractComponent)
    get_param_names(comps::Vector{<:AbstractComponent})

Get the names of parameters from a component or vector of components.
Returns a vector of symbols representing parameter names.
"""
get_param_names(comp::AbstractComponent) = comp.meta.params
get_param_names(comps::Vector{<:AbstractComponent}) = reduce(union, get_param_names.(comps))

"""
    get_nn_names(comp::AbstractComponent)
    get_nn_names(comps::Vector{<:AbstractComponent})

Get the names of neural network components from a component or vector of components.
Returns a vector of symbols representing neural network component names.
"""
get_nn_names(comp::AbstractComponent) = comp.meta.nns
get_nn_names(comps::Vector{<:AbstractComponent}) = reduce(union, get_nn_names.(comps))

"""
    get_var_names(comps::AbstractComponent)

Get all variable names (inputs, outputs, and states) from a component.
Returns a tuple of three vectors: (input_names, output_names, state_names).
"""
get_var_names(comps::AbstractComponent) = get_input_names(comps), get_output_names(comps), get_state_names(comps)

"""
    get_var_names(funcs::Vector{<:AbstractHydroFlux})

Get all variable names from a vector of flux functions, handling dependencies between inputs and outputs.
Returns a tuple of three vectors: (input_names, output_names, state_names).
"""
function get_var_names(funcs::Vector{<:AbstractHydroFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        tmp_input_names = setdiff(get_input_names(func), output_names)
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names
end

"""
    get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractStateFlux})

Get all variable names from vectors of flux and state flux functions, handling dependencies between inputs, outputs, and states.
Returns a tuple of three vectors: (input_names, output_names, state_names).
"""
function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = get_state_names(dfuncs)
    for func in vcat(funcs, dfuncs)
        tmp_input_names = setdiff(setdiff(get_input_names(func), output_names), state_names)
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names, state_names
end

"""
    get_var_names(funcs::Vector{<:AbstractComponent})

Get all variable names from a vector of flux functions, handling dependencies between inputs and outputs.
Returns a tuple of three vectors: (input_names, output_names, state_names).
"""
function get_var_names(components::Vector{<:AbstractComponent})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for comp in components
        tmp_input_names = get_input_names(comp)
        tmp_output_names = get_output_names(comp)
        tmp_state_names = get_state_names(comp)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end
