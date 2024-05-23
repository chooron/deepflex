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
    input::StructArray,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver()
) where {mtk}
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]

    unit_input_names, unit_output_names, unit_state_names = get_var_names(unit)
    unit_var_names = vcat(unit_input_names, unit_output_names, unit_state_names)
    fluxes = StructArray(NamedTuple{Tuple(unit_var_names)}(
        [nm in unit_input_names ? getproperty(input, nm) : zeros(eltype(input[1]), length(timeidx)) for nm in unit_var_names]
    ))
    for tmp_ele in unit.elements
        tmp_output = tmp_ele(fluxes, pas, timeidx=timeidx, solver=solver)
        for nm in vcat(get_var_names(tmp_ele)[[2, 3]]...)
            getproperty(fluxes, nm) .= getproperty(tmp_output, nm)
        end
    end
    fluxes
end

function (unit::HydroUnit{true,false})(
    input::StructArray,
    pas::ComponentVector;
    timeidx::Vector=collect(1:length(input)),
    solver::AbstractSolver=ODESolver()
)
    #* 处理输入pas的异常情况
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]

    #* prepare problem building
    unit_input_names, unit_output_names, unit_state_names = get_var_names(unit)
    unit_var_names = vcat(unit_input_names, unit_output_names, unit_state_names)
    fluxes = StructArray(NamedTuple{Tuple(unit_var_names)}(
        [nm in unit_input_names ? getproperty(input, nm) : zeros(eltype(input[1]), length(timeidx)) for nm in unit_var_names]
    ))
    for nm in unit_state_names
        getproperty(fluxes, nm)[1] = init_states[nm]
    end
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
    fluxes_0 = fluxes[1]
    for tmp_ele in unit.elements
        for tmp_func in tmp_ele.funcs
            tmp_output = NamedTuple{Tuple(get_output_names(tmp_func))}(tmp_func(fluxes_0, params))
            fluxes_0 = merge(fluxes_0, tmp_output)
        end
    end

    # 获取u0和p实际值
    u0 = get_mtk_initstates(elements_systems, fluxes_0, elements_state_names, nfunc_ntp_list)
    p = get_mtk_params(elements_systems, params)
    new_prob = remake(prob, p=p, u0=u0)
    solved_states = solver(new_prob)

    # 更新fluxes中state计算结果
    for (idx, nm) in enumerate(unit_state_names)
        getproperty(fluxes, nm) .= solved_states[idx, :]
    end

    for tmp_ele in unit.elements
        tmp_output = tmp_ele(fluxes, pas, timeidx=timeidx, solved=true)
        for nm in vcat(get_var_names(tmp_ele)[[2, 3]]...)
            getproperty(fluxes, nm) .= getproperty(tmp_output, nm)
        end
    end
    fluxes
end

function (unit::HydroUnit{false,false})(
    input::StructArray,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver()
)
    params = pas[:params] isa Vector ? ComponentVector() : pas[:params]
    init_states = pas[:initstates] isa Vector ? ComponentVector() : pas[:initstates]

    unit_input_names, unit_output_names, unit_state_names = get_var_names(unit)
    unit_var_names = vcat(unit_input_names, unit_output_names, unit_state_names)
    fluxes = StructArray(NamedTuple{Tuple(unit_var_names)}(
        [nm in unit_input_names ? getproperty(input, nm) : zeros(eltype(input[1]), length(timeidx)) for nm in unit_var_names]
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