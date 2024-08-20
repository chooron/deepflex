get_input_names(func::AbstractFlux) = func.infos[:input]
get_input_names(ele::AbstractHydroBucket) = ele.infos[:input]
get_input_names(unit::AbstractModel) = unit.infos[:input]

get_output_names(func::AbstractFlux) = func.infos[:output]
get_output_names(::AbstractStateFlux) = Symbol[]
get_output_names(ele::AbstractHydroBucket) = ele.infos[:output]

get_state_names(::AbstractFlux) = Symbol[]
get_state_names(func::AbstractStateFlux) = [func.infos[:state]]
get_state_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_state_names.(funcs))
get_state_names(ele::AbstractHydroBucket) = ele.infos[:state]

get_param_names(func::AbstractFlux) = func.infos[:param]
get_param_names(::AbstractNeuralFlux) = Symbol[]
get_param_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_param_names.(funcs))
get_param_names(ele::AbstractHydroBucket) = ele.infos[:param]

get_nn_names(::AbstractFlux) = Symbol[]
get_nn_names(func::AbstractNeuralFlux) = func.infos[:param]
get_nn_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_nn_names.(funcs))
get_nn_names(ele::AbstractHydroBucket) = ele.infos[:nn]

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
#* elements name utils
get_var_names(ele::AbstractHydroBucket) = reduce(vcat, collect(ele.infos[[:input, :output, :state]]))
function get_var_names(elements::Vector{<:AbstractHydroBucket})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for ele in elements
        tmp_input_names = get_input_names(ele)
        tmp_output_names = get_output_names(ele)
        tmp_state_names = get_state_names(ele)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end
get_var_names(unit::AbstractModel) = unit.infos[:var]