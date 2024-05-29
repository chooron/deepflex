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
struct HydroUnit{mtk,step} <: AbstractUnit
    "hydrological computation unit name"
    name::Symbol
    "hydrological computation elements"
    elements::Vector{<:AbstractElement}
    "The system of `ModelingToolkit.jl` generated based on multiple elements"
    system::Union{Nothing,ODESystem}

    function HydroUnit(name;
        elements::Vector{<:AbstractElement},
        step::Bool=true)

        mtk = !(false in [typeof(ele).parameters[1] for ele in elements])
        sorted_elements = sort_elements_by_topograph(elements)

        if (!step & mtk)
            system = build_unit_system(sorted_elements, name=name)
        else
            system = nothing
        end

        new{mtk,step}(
            name,
            sorted_elements,
            system,
        )
    end
end

function update_unit!(unit::HydroUnit)
    unit.elements = sort_elements_by_topograph(unit.elements)
end

function add_elements!(unit::HydroUnit; elements::Vector{<:AbstractElement})
    for ele in elements
        push!(unit.elements, ele)
    end
    update_unit!(unit)
end

# 求解并计算
function (unit::HydroUnit{mtk,true})(
    input::NamedTuple,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector,
    solver::AbstractSolver=ODESolver()
) where {mtk}
    fluxes = input
    for tmp_ele in unit.elements
        fluxes = merge(fluxes, tmp_ele(fluxes, pas, timeidx=timeidx, solver=solver))
    end
    fluxes
end

function (unit::HydroUnit{true,false})(
    input::NamedTuple,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector=collect(1:length(input)),
    solver::AbstractSolver=ODESolver()
)
    #* 处理输入pas的异常情况
    params = NamedTuple(pas[:params])
    init_states = NamedTuple(pas[:initstates])

    #* prepare problem building
    unit_input_names, _, unit_state_names = get_var_names(unit)

    nfunc_ntp_list, elements_input_names, elements_state_names = [], [], []
    for ele in unit.elements
        ele_input_names, _, ele_state_names = get_var_names(ele)
        push!(nfunc_ntp_list, extract_neuralflux_ntp(ele.funcs))
        push!(elements_input_names, ele_input_names)
        push!(elements_state_names, ele_state_names)
    end
    elements_systems = [getproperty(unit.system, ele.system.name) for ele in unit.elements]
    #* 问题初始化
    build_sys = setup_input(unit.system, elements_systems, input, timeidx, elements_input_names, unit_input_names, unit.name)
    prob = init_prob(build_sys, elements_systems, nfunc_ntp_list, timeidx)

    #* setup problem with input parameters
    fluxes_0 = merge(NamedTuple{keys(input)}([input[nm][1] for nm in keys(input)]), init_states)
    for tmp_ele in unit.elements
        for tmp_func in tmp_ele.funcs
            tmp_input = [fluxes_0[nm] for nm in get_input_names(tmp_func)]
            tmp_params = [params[nm] for nm in get_param_names(tmp_func)]
            tmp_ouput = NamedTuple{Tuple(get_output_names(tmp_func))}(tmp_func(tmp_input, tmp_params))
            fluxes_0 = merge(fluxes_0, tmp_ouput)
        end
    end
    #* 获取u0和p实际值
    u0 = get_mtk_initstates(elements_systems, fluxes_0, elements_state_names, nfunc_ntp_list)
    p = get_mtk_params(elements_systems, params)
    new_prob = remake(prob, p=p, u0=u0)
    solved_states = solver(new_prob)
    states_ntp = NamedTuple{Tuple(unit_state_names)}(solved_states)

    #* 根据推求的结果计算出其他变量
    fluxes = merge(input, states_ntp)
    for tmp_ele in unit.elements
        fluxes = merge(fluxes, tmp_ele(fluxes, pas, timeidx=timeidx, solved=true))
    end
    fluxes
end

function (unit::HydroUnit{false,false})(
    input::StructArray,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector=collect(1:length(input)),
    solver::AbstractSolver=ODESolver()
)
    #! 当前的求解方式好像无法满足将stateflux合并到一起的求解方式，
    #! 因为后续的stateflux可能需要其他输入变量但是是中间状态
    #* 我想到一个方法或许可以解决这个问题：在构建unit之前就将elements的dflux进行重建
    params = NamedTuple(pas[:params])
    init_states = NamedTuple(pas[:initstates])

    unit_input_names, unit_output_names, unit_state_names = get_var_names(unit)
    unit_var_names = vcat(unit_input_names, unit_output_names, unit_state_names)

    fluxes = StructArray(NamedTuple{Tuple(unit_var_names)}(
        [nm in unit_input_names ? getproperty(input, nm) : zeros(eltype(pas), length(timeidx)) for nm in unit_var_names]
    ))
    for nm in unit_state_names
        @inbounds getproperty(fluxes, nm)[1] = init_states[nm]
    end

    for i in timeidx
        for tmp_ele in unit.elements
            for tmp_func in tmp_ele.funcs
                tmp_output = tmp_func(fluxes[i], params)
                for (idx, nm) in enumerate(get_output_names(tmp_func))
                    @inbounds getproperty(fluxes, nm)[i] = tmp_output[idx]
                end
            end
        end
        if i < length(timeidx)
            for ele in unit.elements
                for dfunc in ele.dfuncs
                    @inbounds getproperty(fluxes, dfunc.state_names)[i+1] = getproperty(fluxes, dfunc.state_names)[i] + dfunc(fluxes[i], params)
                end
            end
        end
    end
    fluxes
end