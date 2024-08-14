get_input_names(func::AbstractFlux) = keys(func.input_info)

get_output_names(func::AbstractFlux) = keys(func.output_info)

get_param_names(func::AbstractFlux) = keys(func.param_info)

get_param_names(::AbstractNeuralFlux) = ()

get_nn_names(::AbstractFlux) = ()

get_nn_names(func::AbstractNeuralFlux) = (func.nn_info[:name],)

function get_input_output_names(funcs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        # 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(get_output_names(func), input_names)
        # 输入需要排除已有的输出变量，表明这个输入是中间计算得到的
        tmp_input_names = setdiff(get_input_names(func), output_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names
end

function get_input_state_names(dfuncs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    state_names = Vector{Symbol}()

    for dfunc in dfuncs
        union!(input_names, get_input_names(dfunc))
        union!(state_names, get_output_names(dfunc))
    end
    input_names, state_names
end

function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux})
    input_names1, output_names = get_input_output_names(funcs)
    input_names2, state_names = get_input_state_names(dfuncs)
    # 避免一些中间变量混淆为输入要素
    setdiff!(input_names2, output_names)
    # 合并两种func的输入要素
    input_names = union(input_names1, input_names2)
    setdiff!(input_names, state_names)
    input_names, output_names, state_names
end

function get_param_names(funcs::Vector{<:AbstractFlux})
    param_names = Vector{Symbol}()
    for func in funcs
        union!(param_names, get_param_names(func))
    end
    param_names
end

get_input_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux}) = get_var_names(funcs, dfuncs)[1]

get_input_names(funcs::Vector{<:AbstractFlux}) = get_input_output_names(funcs)[1]

get_output_names(funcs::Vector{<:AbstractFlux}) = get_input_output_names(funcs)[2]

get_state_names(dfuncs::Vector{<:AbstractFlux}) = get_input_state_names(dfuncs)[2]

get_nn_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_nn_names.(funcs))

#* elements name utils
get_var_names(ele::AbstractHydroElement) = get_var_names(ele.funcs, ele.dfuncs)

get_input_names(ele::AbstractHydroElement) = get_var_names(ele)[1]

get_output_names(ele::AbstractHydroElement) = get_var_names(ele)[2]

get_state_names(ele::AbstractHydroElement) = get_state_names(ele.dfuncs)

get_param_names(ele::AbstractHydroElement) = get_param_names(vcat(ele.funcs, ele.dfuncs))

#* lag element name utils
get_input_names(ele::AbstractLagElement) = collect(reduce(union, get_input_names.(ele.lfuncs)))

get_output_names(ele::AbstractLagElement) = collect(reduce(union, get_output_names.(ele.lfuncs)))

get_param_names(ele::AbstractLagElement) = collect(reduce(union, get_param_names.(ele.lfuncs)))

get_state_names(::AbstractLagElement) = Symbol[]

get_var_names(ele::AbstractLagElement) = get_input_names(ele), get_output_names(ele), get_state_names(ele)

function get_var_names(elements::Vector{<:AbstractElement})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for ele in elements
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(ele)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end

get_input_names(elements::Vector{<:AbstractElement}) = get_var_names(elements)[1]

get_output_names(elements::Vector{<:AbstractElement}) = get_var_names(elements)[2]

get_state_names(elements::Vector{<:AbstractElement}) = reduce(union, [get_state_names(ele) for ele in elements])

get_param_names(elements::Vector{<:AbstractElement}) = reduce(union, [get_param_names(ele) for ele in elements])