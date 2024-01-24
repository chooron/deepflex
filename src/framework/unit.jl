# @kwdef mutable struct Unit{T} <: Component where {T<:Number}
#     id::String

#     # model structure
#     stucture::AbstractGraph

#     # inner variables
#     fluxes::Dict{Symbol,Vector{T}} = Dict()

#     # attribute
#     input_names::Vector{Symbol} = []
# end
function build_unit(; unit::Unit, paraminfos::Vector{ParamInfo})
    parameters = Dict(p.name => p.value for p in paraminfos)
end

function update_fluxes!(unit::Unit; fluxes::ComponentVector)
    unit.fluxes = ComponentVector(unit.fluxes; fluxes...)
end

function get_fluxes(unit::Unit; flux_names::Vector{Symbol})
    output = Dict{Symbol,Vector}()
    for flux_nm in flux_names
        output[flux_nm] = unit.fluxes[flux_nm]
    end
    return ComponentVector(; output...)
end

function get_init_states(sort_eles::Vector{E}) where {E<:BaseElement}
    u_init_dict = Dict{Symbol,Number}()
    for tmp_ele in sort_eles
        if isa(tmp_ele, ODEsElement)
            for sn in tmp_ele.state_names
                u_init_dict[sn] = getproperty(tmp_ele, sn)
            end
        end
    end
    return ComponentVector(u_init_dict)
end


function set_parameters!(unit::Unit; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for idx in topological_sort(unit.structure)
        tmp_ele = get_prop(unit.structure, idx, :ele)
        if isa(tmp_ele, ParameterizedElement) | isa(tmp_ele, StateParameterizedElement)
            set_parameters!(tmp_ele, paraminfos=paraminfos)
        end
        if isa(tmp_ele, StateElement) | isa(tmp_ele, StateParameterizedElement)
            set_states!(tmp_ele, paraminfos=paraminfos)
        end
    end
end

function get_output(unit::Unit; input::ComponentVector{T}, step::Bool=true) where {T<:Number}
    unit_sort_idx = topological_sort(unit.structure)
    unit_sort_eles = [get_prop(unit.structure, idx, :ele) for idx in unit_sort_idx]
    # 开始计算
    if step
        # * This function is calculated element by element
        # initialize unit fluxes
        unit.fluxes = input
        # traversal of the directed graph
        for tmp_ele in unit_sort_eles
            tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
            update_fluxes!(unit, fluxes=tmp_fluxes)
        end
    else
        # * This function is calculated based on the whole Unit
        dt = 1
        xs = 1:dt:length(input[first(keys(input))])
        # fit interpolation functions
        itp = Dict(k => linear_interpolation(xs, input[k]) for k in keys(input))
        # 判断是否存在LuxElement,如果有就是NeuralODE问题
        node = false
        # 获取ODEsElement的所有state初始值
        u_init = get_init_states(unit_sort_eles)
        # 定义整体的ode函数
        function ode_function!(du, u, p, t)
            # element input
            # 使用插值方法获取该时段下的输入值
            tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
            # 遍历Unit中所有的Element进行求解
            for tmp_ele in unit_sort_eles
                # 判断是否为ODEsElement
                if isa(tmp_ele, ODEsElement)
                    # 求解du并更新du
                    tmp_du = get_du(tmp_ele, S=u, input=tmp_input)
                    du = ComponentVector(du; tmp_du...)
                    # 计算出各个flux，更新至tmp_input中
                    tmp_fluxes = get_fluxes(tmp_ele, S=u, input=tmp_input)
                else
                    # 计算出各个flux，更新至tmp_input中
                    tmp_fluxes = get_fluxes(tmp_ele, input=tmp_input)
                end
                tmp_input = ComponentVector(tmp_input; tmp_fluxes...)
            end
        end

        # *求解这个函数
        prob = ODEProblem(ode_function!, u_init, (xs[1], maximum(xs)))
        if node
            sol = solve(prob, BS3(), dt=1.0, saveat=xs, reltol=1e-3, abstol=1e-3, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
        else
            sol = solve(prob, BS3(), p=(), saveat=xs, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
        end
        solved_u = sol.u
        solved_u_matrix = hcat(solved_u...)
        solved_u = ComponentVector(; Dict(nm => solved_u_matrix[idx, :] for (idx, nm) in enumerate(keys(solved_u[1])))...)
        unit.fluxes = input
        # 带入求解的结果计算最终的输出结果
        for tmp_ele in unit_sort_eles
            if isa(tmp_ele, ODEsElement)
                tmp_fluxes = get_fluxes(tmp_ele, S=solved_u, input=unit.fluxes)
            else
                tmp_fluxes = get_fluxes(tmp_ele, input=unit.fluxes)
            end
            update_fluxes!(unit, fluxes=tmp_fluxes)
        end
    end
    return get_fluxes(unit, flux_names=unit.output_names)
end
