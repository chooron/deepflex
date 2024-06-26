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
struct HydroUnit{mtk} <: AbstractUnit
    "hydrological computation unit name"
    name::Symbol
    "hydrological computation elements"
    elements::Vector{<:AbstractElement}

    function HydroUnit(name; elements::Vector{<:AbstractElement})

        mtk = !(false in [typeof(ele).parameters[1] for ele in elements])
        # sorted_elements = sort_elements_by_topograph(elements)
        # println([ele.name for ele in sorted_elements])

        # if mtk
        #     system = build_unit_system(sorted_elements, name=name)
        # else
        #     system = nothing
        # end

        new{mtk}(
            name,
            elements,
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
function (unit::HydroUnit)(
    input::NamedTuple,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector,
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    for tmp_ele in unit.elements
        fluxes = merge(fluxes, tmp_ele(fluxes, pas, timeidx=timeidx, solver=solver))
    end
    fluxes
end
