#* name utils for flux
function get_input_names(func::Union{AbstractSimpleFlux,AbstractNeuralFlux})
    if eltype(func.input_names) isa Pair
        input_names = [v for (_, v) in func.input_names]
    elseif func.input_names isa Vector
        input_names = func.input_names
    else
        input_names = [func.input_names]
    end
    input_names
end

function get_output_names(func::Union{AbstractSimpleFlux,AbstractNeuralFlux})
    if func.output_names isa Vector
        output_names = func.output_names
    else
        output_names = [func.output_names]
    end
    output_names
end

function get_param_names(func::AbstractFlux)
    if func.param_names isa Vector
        param_names = func.param_names
    else
        param_names = [func.param_names]
    end
    param_names
end

function get_input_names(func::AbstractStateFlux)
    vcat(func.influx_names, func.outflux_names)
end

function get_output_names(func::AbstractStateFlux)
    func.state_names
end

function get_input_names(func::AbstractLagFlux)
    func.flux_name
end

function get_output_names(func::AbstractLagFlux)
    func.flux_name
end

function get_param_names(func::AbstractLagFlux)
    func.lag_time
end

function get_func_io_names(funcs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        # extract the input and output name in flux
        # 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(get_output_names(func), input_names)
        # 输入需要排除已有的输出变量，表明这个输入是中间计算得到的
        tmp_input_names = setdiff(get_input_names(func), output_names)
        # 合并名称
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names
end

function get_dfunc_io_names(dfuncs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    state_names = Vector{Symbol}()

    for dfunc in dfuncs
        union!(input_names, get_input_names(dfunc))
        union!(state_names, get_output_names(dfunc))
    end
    input_names, state_names
end

function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux})
    input_names1, output_names = get_func_io_names(funcs)
    input_names2, state_names = get_dfunc_io_names(dfuncs)
    # 避免一些中间变量混淆为输入要素
    setdiff!(input_names2, output_names)
    # 合并两种func的输入要素
    input_names = union(input_names1, input_names2)
    setdiff!(input_names, state_names)
    input_names, output_names, state_names
end

get_input_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux}) = get_var_names(funcs, dfuncs)[1]

get_input_names(funcs::Vector{<:AbstractFlux}) = get_func_io_names(funcs)[1]

get_output_names(funcs::Vector{<:AbstractFlux}) = get_func_io_names(funcs)[2]

get_state_names(dfuncs::Vector{<:AbstractFlux}) = get_dfunc_io_names(dfuncs)[2]


function get_param_names(funcs::Vector{<:AbstractFlux})
    param_names = Vector{Symbol}()
    for func in funcs
        union!(param_names, get_param_names(func))
    end
    param_names
end

#* name utils for elements
function get_ele_io_names(elements::Vector{<:AbstractElement})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for ele in elements
        union!(input_names, setdiff(get_input_names(ele.funcs, ele.dfuncs), output_names))
        union!(output_names, get_output_names(ele.funcs))
    end
    input_names, output_names
end

function get_ele_param_names(elements::Vector{<:AbstractElement})
    param_names = Vector{Symbol}()
    for ele in elements
        union!(param_names, get_param_names(vcat(ele.funcs, ele.dfuncs)))
    end
    param_names
end

function get_ele_state_names(elements::Vector{<:AbstractElement})
    state_names = Vector{Symbol}()
    for ele in elements
        union!(state_names, get_state_names(ele.dfuncs))
    end
    state_names
end

#* name utils for elements
function get_ele_io_names(units::NamedTuple)
    input_names_list = []
    output_names_list = []
    for unit in keys(units)
        ele_input_names, ele_output_names = get_ele_io_names(unit)
        push!(input_names_list, ele_input_names)
        push!(output_names_list, ele_output_names)
    end
    namedtuple(keys(units), input_names_list), namedtuple(keys(units), output_names_list)
end

function get_ele_param_names(units::NamedTuple)
    param_tuple_list = []
    for nm in keys(units)
        push!(param_tuple_list, (unit=get_param_names(units[nm]), route=get_param_names([route])))
    end
    param_tuple_list
end

function get_ele_state_names(units::NamedTuple)
    state_tuple_list = []
    for nm in keys(units)
        push!(state_tuple_list, (unit=get_state_names(units[nm]), route=get_state_names([route])))
    end
    state_tuple_list
end