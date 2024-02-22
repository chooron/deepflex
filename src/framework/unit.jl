@kwdef mutable struct Unit{T,E} <: AbstractUnit where {T<:Number,E<:AbstractElement}
    name::String

    # parameters
    parameters::ComponentVector{T}

    # model structure
    elements::Vector{E}

    # attribute
    input_names::Set{Symbol}
    output_names::Set{Symbol}
    state_names::Set{Symbol}

    # inner variables
    fluxes::ComponentVector = ComponentVector()
end

function build_unit(; name::String, elements::Vector{E}) where {E<:AbstractElement}
    parameters = ComponentVector()
    input_names = Set{Symbol}()
    output_names = Set{Symbol}()
    state_names = Set{Symbol}()
    for ele in elements
        parameters = ComponentVector(parameters; ele.parameters...)
        union!(input_names, ele.input_names)
        union!(output_names, ele.output_names)
        union!(state_names, ele.state_names)
    end
    Unit(
        name=name,
        parameters=parameters,
        elements=elements,
        input_names=input_names,
        output_names=output_names,
        state_names=state_names,
    )
end

function pretrain!(unit::AbstractUnit; input::ComponentVector{T}, train_config...) where {T<:Number}
    for ele in unit.elements
        pretrain!(ele; input=input, train_config...)
    end
end

function get_fluxes(unit::AbstractUnit; flux_names::Set{Symbol})
    output = Dict{Symbol,Vector}()
    for flux_nm in flux_names
        output[flux_nm] = unit.fluxes[flux_nm]
    end
    return ComponentVector(; output...)
end

function get_states(unit::AbstractUnit; state_names::Union{Set{Symbol},Nothing}=nothing)
    states = ComponentVector()
    for ele in unit.elements
        states = ComponentVector(states; get_states(ele, state_names=state_names)...)
    end
    return states
end

function get_init_states(sort_eles::Vector{E}) where {E<:AbstractElement}
    init_states = ComponentVector()
    for ele in sort_eles
        init_states = ComponentVector(init_states; ele.init_states...)
    end
    return init_states
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

function get_output(unit::Unit; input::ComponentVector{T}, step::Bool=true, kwargs...) where {T<:Number}
    # 开始计算
    if step
        # * This function is calculated element by element
        # initialize unit fluxes
        unit.fluxes = input
        # traversal of the directed graph
        for tmp_ele in unit.elements
            tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
            unit.fluxes = ComponentVector(unit.fluxes; tmp_fluxes...)
        end
    else
        # * This function is calculated based on the whole Unit
        dt = 1
        xs = 1:dt:length(input[first(keys(input))])
        # fit interpolation functions
        itp = Dict(k => linear_interpolation(xs, input[k]) for k in keys(input))
        # 获取ODEsElement的所有state初始值
        u_init = get_init_states(unit.elements)
        # 定义整体的ode函数
        function ode_function!(du, u, p, t)
            # element input
            # 使用插值方法获取该时段下的输入值
            tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
            # 遍历Unit中所有的Element进行求解
            for tmp_ele in unit.elements
                # 计算出各个flux，更新至tmp_input中
                tmp_fluxes = get_fluxes(tmp_ele, state=u, input=tmp_input)
                # 求解du并更新du
                tmp_du = tmp_ele.get_du(tmp_fluxes, tmp_ele.parameters)
                for k in tmp_ele.state_names
                    du[k] = tmp_du[k]
                end
                tmp_input = ComponentVector(tmp_input; tmp_fluxes...)
            end
        end

        # *求解这个函数
        prob = ODEProblem(ode_function!, u_init, (xs[1], maximum(xs)))
        sol = solve(
            prob,
            get(kwargs, :alg, BS3()),
            saveat=xs,
            dt=dt,
            reltol=get(kwargs, :reltol, 1e-3),
            abstol=get(kwargs, :abstol, 1e-3),
            sensealg=get(kwargs, :sensealg, ForwardDiffSensitivity())
        )
        solved_u = sol.u
        solved_u_matrix = hcat(solved_u...)
        solved_u = ComponentVector(; Dict(nm => solved_u_matrix[idx, :] for (idx, nm) in enumerate(keys(solved_u[1])))...)
        for ele in unit.elements
            ele.states = solved_u[collect(ele.state_names)]
        end
        unit.fluxes = input
        # *带入求解的结果计算最终的输出结果
        for tmp_ele in unit.elements
            tmp_fluxes = get_fluxes(tmp_ele, state=solved_u, input=unit.fluxes)
            unit.fluxes = ComponentVector(unit.fluxes; tmp_fluxes...)
        end
    end
    # return get_fluxes(unit, flux_names=unit.output_names)
    return unit.fluxes
end
