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
struct HydroUnit <: AbstractUnit
    "hydrological computation unit name"
    name::Symbol
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}

    function HydroUnit(name; components::Vector{<:AbstractComponent})
        new(
            name,
            components,
        )
    end
end

# 求解并计算
function (unit::HydroUnit)(
    input::NamedTuple,
    pas::Union{ComponentVector,NamedTuple};
    timeidx::Vector,
    solver::AbstractSolver=ODESolver()
)
    fluxes = input
    for tmp_ele in unit.components
        tmp_fluxes = tmp_ele(fluxes, pas, timeidx=timeidx, solver=solver)
        fluxes = merge(fluxes, tmp_fluxes)
    end
    fluxes
end
