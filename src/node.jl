struct HydroNode{step} <: AbstractComponent
    name::Symbol
    #* 子单元名称
    subnames::Vector{Symbol}
    #* 垂向计算层
    layers::Vector{Vector{<:AbstractElement}}
    #* 横向计算层(坡地产流模块)
    routes::Vector{<:AbstractElement}
    #* 节点对应的面积
    area::Number
    #* 根据mtk生成的system
    sys::Union{Nothing,Vector{ODESystem}}

    function HydroNode(name; layers::NamedTuple, routes::NamedTuple, area::Number=100.0, step::Bool=true)
        unit_names = keys(layers)
        if step
            sys_list = nothing
        else
            sys_list = [build_node_system(layers[nm], name=nm) for nm in unit_names]
        end
        new{step}(
            name,
            collect(keys(layers)),
            collect(layers),
            collect(routes),
            area,
            sys_list
        )
    end
end

# todo parallel computing
function (node::HydroNode{true})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for (idx, nm) in enumerate(node.subnames)
        tmp_input = input[nm]
        tmp_unit_ouput = node.layers[idx](tmp_input, pas[nm], solver=solver)
        tmp_route_output = node.routes[idx](tmp_unit_ouput, pas[nm])
        push!(output_list, tmp_route_output[:flow] .* pas[nm][:weight])
    end
    (flow=sum(output_list),)
end

function (node::HydroNode{false})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    output_list = []
    for (idx, nm) in enumerate(node.subnames)
        tmp_input = input[nm]
        prob = setup_input(node.layers[idx], node.sys[idx], input=tmp_input, name=nm)
        new_prob = setup_prob(node.layers[idx], node.sys[idx], prob, params=pas[nm][:params], init_states=pas[nm][:initstates])
        solved_states = solver(new_prob, get_ele_state_names(node.layers[idx]))
        tmp_input = merge(tmp_input, solved_states)
        tmp_unit_ouput = node.layers[idx](tmp_input, pas[nm], solver=solver)
        tmp_route_output = node.routes[idx](tmp_unit_ouput, pas[nm])
        push!(output_list, tmp_route_output[:flow] .* pas[nm][:weight])
    end
    (flow=sum(output_list),)
end

function build_node_system(
    elements::AbstractVector{HydroElement{true}};
    name::Symbol,
)
    eqs = []
    elements = elements
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in elements
        #* 连接element之间的变量
        for tmp_ele2 in filter(ele -> ele.name != tmp_ele1.name, elements)
            share_var_names = intersect(
                vcat(get_var_names(tmp_ele1.funcs, tmp_ele1.dfuncs)...),
                vcat(get_var_names(tmp_ele2.funcs, tmp_ele2.dfuncs)...)
            )
            for nm in share_var_names
                push!(eqs, getproperty(tmp_ele1.sys, nm) ~ getproperty(tmp_ele2.sys, nm))
            end
        end
    end
    compose(ODESystem(eqs, t; name=Symbol(name, :_sys)), [ele.sys for ele in elements]...)
end

function setup_input(
    layers::AbstractVector{HydroElement{true}},
    sys::ODESystem;
    input::NamedTuple,
    name::Symbol,
)
    #* 首先构建data的插值系统
    eqs = Equation[]
    node_varinfo = merge([ele.varinfo for ele in layers]...)
    unit_input_names = get_ele_io_names(layers)[1]
    itp_sys = build_itp_system(NamedTupleTools.select(input, unit_input_names), input[:time], node_varinfo, name=name)
    for ele in layers
        for nm in filter(nm -> nm in get_input_names(ele.funcs, ele.dfuncs), keys(input))
            push!(eqs, getproperty(getproperty(sys, ele.sys.name), nm) ~ getproperty(itp_sys, nm))
        end
    end
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(name, :_comp_sys)), sys, itp_sys)
    sys = structural_simplify(compose_sys)
    prob = ODEProblem(sys, Pair[], (input[:time][1], input[:time][end]), [], warn_initialize_determined=true)
    prob
end

function setup_prob(
    layers::AbstractVector{HydroElement{true}},
    sys::ODESystem,
    prob::ODEProblem;
    params::ComponentVector,
    init_states::ComponentVector,
)
    #* setup init states
    u0 = Pair[]
    #* setup parameters
    p = Pair[]
    for ele in layers
        tmp_ele_sys = getproperty(sys, ele.sys.name)
        for nm in filter(nm -> nm in get_state_names(ele.dfuncs), keys(init_states))
            push!(u0, getproperty(tmp_ele_sys, Symbol(nm)) => init_states[Symbol(nm)])
        end
        for nm in ModelingToolkit.parameters(ele.sys)
            push!(p, getproperty(tmp_ele_sys, Symbol(nm)) => params[Symbol(nm)])
        end
    end
    new_prob = remake(prob, p=p, u0=u0)
    new_prob
end