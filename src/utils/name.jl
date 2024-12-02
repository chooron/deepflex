
"""
Metadata about the component, including:
- name: Symbol representing the component's name
- inputs: Vector of input variable names
- outputs: Vector of output variable names
- params: Vector of parameter names
- states: Vector of state variable names
- nns: Vector of neural network names
"""
struct HydroMeta
    name::Symbol
    inputs::Vector{Symbol}
    outputs::Vector{Symbol}
    params::Vector{Symbol}
    states::Vector{Symbol}
    nns::Vector{Symbol}

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

get_input_names(func::AbstractFlux) = func.meta.inputs
get_input_names(ele::AbstractBucket) = ele.meta.inputs
get_input_names(route::AbstractRoute) = route.meta.inputs
get_input_names(unit::AbstractModel) = unit.meta.inputs
get_input_names(comp::AbstractComponent) = comp.meta.inputs

get_output_names(func::AbstractFlux) = func.meta.outputs
get_output_names(::AbstractStateFlux) = Symbol[]
get_output_names(ele::AbstractBucket) = ele.meta.outputs
get_output_names(route::AbstractRoute) = route.meta.outputs
get_output_names(unit::AbstractModel) = unit.meta.outputs
get_output_names(comp::AbstractComponent) = comp.meta.outputs

get_state_names(flux::AbstractFlux) = flux.meta.states
get_state_names(fluxes::Vector{<:AbstractFlux}) = reduce(union, get_state_names.(fluxes))
get_state_names(ele::AbstractBucket) = ele.meta.states
get_state_names(route::AbstractRoute) = route.meta.states
get_state_names(unit::AbstractModel) = unit.meta.states
get_state_names(comp::AbstractComponent) = comp.meta.states

get_param_names(func::AbstractFlux) = func.meta.params
get_param_names(::AbstractNeuralFlux) = Symbol[]
get_param_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_param_names.(funcs))
get_param_names(ele::AbstractBucket) = ele.meta.params
get_param_names(route::AbstractRoute) = route.meta.params
get_param_names(unit::AbstractModel) = unit.meta.params
get_param_names(comp::AbstractComponent) = comp.meta.params

get_nn_names(::AbstractFlux) = Symbol[]
get_nn_names(func::AbstractNeuralFlux) = func.meta.nns
get_nn_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_nn_names.(funcs))
get_nn_names(ele::AbstractBucket) = ele.meta.nns
get_nn_names(route::AbstractRoute) = route.meta.nns
get_nn_names(unit::AbstractModel) = unit.meta.nns
get_nn_names(comp::AbstractComponent) = comp.meta.nns

get_var_names(func::AbstractFlux) = get_input_names(func), get_output_names(func), get_state_names(func)
#* elements name utils
get_var_names(ele::AbstractBucket) = ele.meta.inputs, ele.meta.outputs, ele.meta.states
get_var_names(route::AbstractRoute) = route.meta.inputs, route.meta.outputs, route.meta.states
get_var_names(unit::AbstractModel) = unit.meta.inputs, unit.meta.outputs, unit.meta.states

function get_var_names(funcs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        #* 输入需要排除已有的输出变量，表明这个输入是中间计算得到的，此外state也不作为输入
        tmp_input_names = setdiff(get_input_names(func), output_names)
        #* 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names
end

function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = get_state_names(dfuncs)
    for func in vcat(funcs, dfuncs)
        #* 输入需要排除已有的输出变量，表明这个输入是中间计算得到的，此外state也不作为输入
        tmp_input_names = setdiff(setdiff(get_input_names(func), output_names), state_names)
        #* 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names, state_names
end

function get_var_names(components::Vector{<:AbstractComponent})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for com in components
        tmp_input_names = get_input_names(com)
        tmp_output_names = get_output_names(com)
        tmp_state_names = get_state_names(com)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end
