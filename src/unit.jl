mutable struct HydroUnit{step} <: AbstractUnit
    #* 单元名称
    name::Symbol
    ##* 一个基础单元必须包括的计算层
    #* 地表水层Surface, 典型Element包括降雪模块，截流模块，蒸发模块，下渗模块等
    surface::Vector{<:AbstractElement}
    #* 壤中水层Soil，典型Element包括土壤水分模块，产流计算模块等
    soil::Vector{<:AbstractElement}
    #* 自由水层FreeWater，典型Element包括地下水，地表水，壤中流等
    freewater::Vector{<:AbstractElement}
    #* 针对多个elements的计算图
    topology::Union{Nothing,SimpleDiGraph}
    ##* 根据mtk生成的system
    system::Union{Nothing,ODESystem}

    function HydroUnit(name;
        surface::Union{AbstractElement,Vector{<:AbstractElement}},
        soil::Union{AbstractElement,Vector{<:AbstractElement}},
        freewater::Union{AbstractElement,Vector{<:AbstractElement}},
        step::Bool=true)
        if step
            topology = build_compute_topology(vcat(surface, soil, freewater))
            sys = nothing
        else
            topology = nothing
            sys = build_unit_system(vcat(surface, soil, freewater), name=name)
        end
        new{step}(
            name,
            surface,
            soil,
            freewater,
            topology,
            sys
        )
    end
end

function update_attr!(unit::HydroUnit{true})
    unit.topology = build_compute_topology(vcat(unit.surface, unit.soil, unit.freewater))
end

function update_attr!(unit::HydroUnit{false})
    unit.system = build_unit_system(vcat(surface, soil, freewater), name=unit.name)
end

function add_elements!(unit::HydroUnit; elements::Vector{<:AbstractElement})
    for ele in elements
        if ele.etype == SurfaceType
            push!(unit.surface, ele)
        elseif ele.etype == SoilType
            push!(unit.soil, ele)
        elseif ele.etype == FreeWaterType
            push!(unit.freewater, ele)
        end
    end
    update_attr!(unit)
end

function remove_elements!(unit::HydroUnit; elements::Vector{Symbol})
    #! to be implemented
end

# todo 需要提供任意层都能够随意添加计算模块的功能，因此这个同样也需要网络计算
function (elements::AbstractVector{<:HydroElement})(
    input::NamedTuple,
    pas::ComponentVector;
    topology::SimpleDiGraph,
    solver::AbstractSolver=ODESolver()
)
    #* 构建elements map信息
    ele_output_and_state_names(ele) = vcat(get_output_names(ele.funcs), get_state_names(ele.dfuncs))
    elements_ntp = namedtuple(
        vcat([ele_output_and_state_names(ele) for ele in elements]...),
        vcat([repeat([ele], length(ele_output_and_state_names(ele))) for ele in elements]...)
    )
    elements_var_names = unique(vcat(get_ele_io_names(elements)..., get_ele_state_names(elements)))

    #* 针对多个element的网络计算
    fluxes = input
    for flux_idx in topological_sort(topology)
        tmp_flux_name = elements_var_names[flux_idx]
        if !(tmp_flux_name in keys(fluxes))
            tmp_ele = elements_ntp[tmp_flux_name]
            if length(tmp_ele.dfuncs) > 0
                solve_states = solve_prob(tmp_ele, input=fluxes, pas=pas, solver=solver)
                fluxes = merge(fluxes, solve_states)
            end
            fluxes = merge(fluxes, tmp_ele(fluxes, pas))
        end
    end
    fluxes
end

# 求解并计算
function (unit::HydroUnit{true})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    for (elements, topology) in zip(get_all_elements(unit), get_all_topologys(unit))
        tmp_fluxes = elements(fluxes, pas, topology=topology, solver=solver)
        fluxes = merge(fluxes, tmp_fluxes)
    end
    fluxes
end

function (unit::HydroUnit{false})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    all_elements = get_all_elements(unit)
    prob = setup_input(all_elements, unit.system, input=fluxes, name=name)
    new_prob = setup_prob(all_elements, unit.system, prob, params=pas[:params], init_states=pas[:initstates])
    solved_states = solver(new_prob, get_ele_state_names(all_elements))
    fluxes = merge(fluxes, solved_states)
    for ele in all_elements
        fluxes = merge(fluxes, ele(fluxes, pas))
    end
    fluxes
end

function build_unit_system(
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
    elements::AbstractVector{HydroElement{true}},
    sys::ODESystem;
    input::NamedTuple,
    name::Symbol,
)
    #* 首先构建data的插值系统
    eqs = Equation[]
    node_varinfo = merge([ele.varinfo for ele in elements]...)
    unit_input_names = get_ele_io_names(elements)[1]
    itp_sys = build_itp_system(NamedTupleTools.select(input, unit_input_names), input[:time], node_varinfo, name=name)
    for ele in elements
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
    elements::AbstractVector{HydroElement{true}},
    sys::ODESystem,
    prob::ODEProblem;
    params::ComponentVector,
    init_states::ComponentVector,
)
    #* setup init states
    u0 = Pair[]
    #* setup parameters
    p = Pair[]
    for ele in elements
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