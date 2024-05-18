"""
$(TYPEDEF)
The basic hydrological calculation unit usually contains multiple hydrological calculation modules,
Each hydrological calculation module can be solved through step-by-step calculation and overall calculation.
It is usually used to simulate the vertical calculation process of a watershed or a calculation unit.

a basic hydrology unit must include:
    Surface water layer, typical Elements include snowfall module, interception module, evaporation module, infiltration module, etc.
    Soil, the water layer in the soil. Typical Elements include soil moisture module, runoff calculation module, etc.
    FreeWater layer, typical elements include groundwater, surface water, soil flow, etc.

# Fields
$(FIELDS)
# Example
```
using LumpedHydro
using LumpedHydro.ExpHydro: Surface, Soil, FreeWater

HydroUnit(
    name,
    elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk), FreeWater(name=name, mtk=mtk)],
    step=step,
)
```
"""
mutable struct HydroUnit{mtk,step} <: AbstractUnit
    "hydrological computation unit name"
    name::Symbol
    "hydrological computation elements"
    elements::Vector{<:AbstractElement}
    "Computational graph for multiple elements"
    topology::MetaTopology
    "The system of `ModelingToolkit.jl` generated based on multiple elements"
    system::Union{Nothing,ODESystem}

    function HydroUnit(name;
        elements::Vector{<:AbstractElement},
        step::Bool=true)

        topology = build_compute_topology(elements)
        mtk = !(false in [typeof(ele).parameters[1] for ele in elements])
        if (!step & mtk)
            system = build_unit_system(elements, name=name)
        else
            system = nothing
        end

        new{mtk,step}(
            name,
            elements,
            topology,
            system,
        )
    end
end

function update_unit!(unit::HydroUnit{true})
    unit.topology = build_compute_topology(unit.elements)
end

function update_unit!(unit::HydroUnit{false})
    unit.system = build_unit_system(unit.elements, name=unit.name)
end

function add_elements!(unit::HydroUnit; elements::Vector{<:AbstractElement})
    for ele in elements
        push!(unit.elements, ele)
    end
    update_attr!(unit)
end

function remove_elements!(unit::HydroUnit; elements::Vector{Symbol})
    #! to be implemented
end

# 求解并计算
function (unit::HydroUnit{mtk,true})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
) where {mtk}
    fluxes = input
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]
    #* 先求解ODEProblem
    for flux_idx in topological_sort(unit.topology.digraph)
        tmp_flux_name = unit.topology.node_names[flux_idx]
        if !(tmp_flux_name in keys(fluxes))
            tmp_ele = unit.topology.node_maps[tmp_flux_name]
            if length(tmp_ele.dfuncs) > 0
                solved_states = solve_prob(tmp_ele, input=fluxes,
                    params=params, init_states=init_states, solver=solver)
                fluxes = merge(fluxes, solved_states)
            end
            fluxes = merge(fluxes, tmp_ele(fluxes, params))
        end
    end
    fluxes
end

function (unit::HydroUnit{false,false})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]
    init_states_names = collect(keys(init_states))

    unit_input_names = get_input_names(unit)
    unit_all_dfuncs = vcat([ele.dfuncs for ele in unit.elements]...)
    itp_dict = namedtuple(unit_input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in unit_input_names])

    function single_unit_ode_func!(du, u, p, t)
        tmp_fluxes = namedtuple(unit_input_names, [itp_dict[nm](t) for nm in unit_input_names])
        tmp_states = namedtuple(init_states_names, u)
        tmp_fluxes = merge(tmp_fluxes, tmp_states)
        tmp_fluxes = unit.elements(tmp_fluxes, params, unit.topology)
        du[:] = [dfunc(tmp_fluxes, params)[dfunc.state_names] for dfunc in unit_all_dfuncs]
    end

    prob = ODEProblem(single_unit_ode_func!, collect(init_states), (input[:time][1], input[:time][end]))
    solved_states = solver(prob, init_states_names)
    fluxes = merge(input, solved_states)
    unit.elements(fluxes, params, unit.topology)
end


function (unit::HydroUnit{true,false})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    #* 处理输入pas的异常情况
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]

    #* prepare problem building
    fluxes = input
    unit_input_names = get_input_names(unit)
    nfunc_ntp_list = [extract_neuralflux_ntp(ele.funcs) for ele in unit.elements]
    elements_input_names = [get_input_names(ele) for ele in unit.elements]
    elements_state_names = [get_state_names(ele) for ele in unit.elements]
    elements_systems = [getproperty(unit.system, ele.system.name) for ele in unit.elements]

    #* 问题初始化
    build_sys = setup_input(
        unit.system,
        elements_systems,
        input,
        elements_input_names,
        unit.name
    )
    prob = init_prob(build_sys, elements_systems, nfunc_ntp_list, input[:time])

    #* setup problem with input parameters
    input_0 = namedtuple(unit_input_names, [input[nm][1] for nm in unit_input_names])
    initstates_ntp = namedtuple(keys(init_states), [init_states[nm] for nm in keys(init_states)])
    sol_0 = unit.elements(merge(input_0, initstates_ntp), params, unit.topology)

    # 获取u0和p实际值
    u0 = get_mtk_initstates(elements_systems, sol_0, elements_state_names, nfunc_ntp_list)
    p = get_mtk_params(elements_systems, params)

    new_prob = remake(prob, p=p, u0=u0)

    solved_states = solver(new_prob, get_state_names(unit.elements))
    fluxes = merge(fluxes, solved_states)
    fluxes = unit.elements(fluxes, params, unit.topology)
    fluxes
end

function build_unit_system(
    elements::AbstractVector{<:HydroElement};
    name::Symbol,
)
    eqs = []
    ele_input_names = get_input_names(elements)
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in elements
        #* 连接element之间的变量
        for tmp_ele2 in filter(ele -> ele.name != tmp_ele1.name, elements)
            #* 这里是为了找到element之间共享的flux但是不是作为unit输入的flux
            share_var_names = intersect(vcat(get_var_names(tmp_ele1)...), vcat(get_var_names(tmp_ele2)...))
            for nm in setdiff(share_var_names, ele_input_names)
                push!(eqs, getproperty(tmp_ele1.system, nm) ~ getproperty(tmp_ele2.system, nm))
            end
        end
    end
    compose(ODESystem(eqs, t; name=Symbol(name, :_sys)), [ele.system for ele in elements]...)
end