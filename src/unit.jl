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
mutable struct HydroUnit{step} <: AbstractUnit
    "hydrological computation unit name"
    name::Symbol
    "hydrological computation elements"
    elements::Vector{<:AbstractElement}
    "Computational graph for multiple elements"
    topology::Union{Nothing,SimpleDiGraph}
    "The system of `ModelingToolkit.jl` generated based on multiple elements"
    system::Union{Nothing,ODESystem}

    function HydroUnit(name;
        elements::Vector{<:AbstractElement},
        step::Bool=true)

        topology = build_compute_topology(elements)

        if step
            system = nothing
        else
            system = build_unit_system(elements, name=name)
        end

        new{step}(
            name,
            elements,
            topology,
            system
        )
    end
end

function update_attr!(unit::HydroUnit{true})
    unit.topology = build_compute_topology(unit.elements)
end

function update_attr!(unit::HydroUnit{false})
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

function (elements::AbstractVector{<:HydroElement})(
    input::NamedTuple,
    pas::ComponentVector;
    topology::SimpleDiGraph,
    solver::AbstractSolver=ODESolver()
)
    #* 构建elements map信息
    ele_output_and_state_names(ele) = vcat(get_output_names(ele), get_state_names(ele))
    elements_ntp = namedtuple(
        vcat([ele_output_and_state_names(ele) for ele in elements]...),
        vcat([repeat([ele], length(ele_output_and_state_names(ele))) for ele in elements]...)
    )
    elements_var_names = unique(vcat(get_input_output_names(elements)..., get_state_names(elements)))

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
    tmp_fluxes = unit.elements(fluxes, pas, topology=unit.topology, solver=solver)
    merge(fluxes, tmp_fluxes)
end

function (unit::HydroUnit{false})(
    input::NamedTuple,
    pas::ComponentVector;
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    prob = setup_input(unit, input=fluxes, name=unit.name)
    new_prob = setup_prob(unit, prob, input=fluxes, params=pas[:params], init_states=pas[:initstates])
    solved_states = solver(new_prob, get_state_names(unit.elements))
    fluxes = merge(fluxes, solved_states)
    fluxes = unit.elements(fluxes, pas, topology=unit.topology, solver=solver)
    fluxes
end

function build_unit_system(
    elements::AbstractVector{<:HydroElement};
    name::Symbol,
)
    eqs = []
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in elements
        #* 连接element之间的变量
        for tmp_ele2 in filter(ele -> ele.name != tmp_ele1.name, elements)
            share_var_names = intersect(vcat(get_var_names(tmp_ele1)...), vcat(get_var_names(tmp_ele2)...))
            for nm in share_var_names
                push!(eqs, getproperty(tmp_ele1.sys, nm) ~ getproperty(tmp_ele2.sys, nm))
            end
        end
    end
    compose(ODESystem(eqs, t; name=Symbol(name, :_sys)), [ele.sys for ele in elements]...)
end

function setup_input(
    unit::HydroUnit;
    input::NamedTuple,
    name::Symbol,
)
    #* 首先构建data的插值系统
    eqs = Equation[]
    elements_varinfo = merge([ele.varinfo for ele in unit.elements]...)
    unit_input_names = get_input_names(unit)
    itp_sys = build_itp_system(NamedTupleTools.select(input, unit_input_names), input[:time], elements_varinfo, name=name)
    for ele in unit.elements
        tmp_ele_sys = getproperty(unit.system, ele.sys.name)
        for nm in filter(nm -> nm in get_input_names(ele), keys(input))
            push!(eqs, getproperty(tmp_ele_sys, nm) ~ getproperty(itp_sys, nm))
        end
    end
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(name, :_comp_sys)), unit.system, itp_sys)
    #* 设置系统初始值并构建problem
    build_u0 = Pair[]
    for ele in unit.elements
        tmp_ele_sys = getproperty(unit.system, ele.sys.name)
        for func in filter(func -> func isa AbstractNeuralFlux, ele.funcs)
            func_nn_sys = getproperty(tmp_ele_sys, func.chain_name)
            push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(get_input_names(func))))
        end
    end
    build_sys = structural_simplify(compose_sys)
    prob = ODEProblem(build_sys, build_u0, (input[:time][1], input[:time][end]), [], warn_initialize_determined=true)
    prob
end

function setup_prob(
    unit::HydroUnit,
    prob::ODEProblem;
    input::NamedTuple,
    params::ComponentVector,
    init_states::ComponentVector,
    kw...
)
    #* setup init states
    u0 = Pair[]
    #* setup parameters
    p = Pair[]
    for ele in unit.elements
        tmp_ele_sys = getproperty(unit.system, ele.sys.name)
        tmp_sol_0 = get_sol_0(ele, input=input, params=params, init_states=init_states)
        nfunc_ntp = extract_neuralflux_ntp(ele.funcs)
        tmp_u0 = get_mtk_initstates(tmp_ele_sys; sol_0=tmp_sol_0, state_names=get_state_names(ele), nfunc_ntp=nfunc_ntp)
        tmp_p = get_mtk_params(tmp_ele_sys; params=params)
        u0 = vcat(u0, tmp_u0)
        p = vcat(p, tmp_p)
    end
    new_prob = remake(prob, p=p, u0=u0)
    new_prob
end