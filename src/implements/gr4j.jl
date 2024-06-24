@reexport module GR4J

using ..LumpedHydro

"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:prcp, :pet] => [:rainfall]),
        SimpleFlux([:rainfall] => [:infiltration])
    ]

    HydroElement(
        Symbol(name, :_surface),
        funcs=funcs,
        mtk=mtk
    )
end

"""
SoilWaterReservoir in GR4J
"""
function Soil(; name::Symbol, mtk::Bool=true)

    fluxes = [
        SimpleFlux([:soilwater, :infiltration] => [:saturation], [:x1]),
        SimpleFlux([:soilwater, :pet] => [:evap], [:x1]),
        SimpleFlux([:soilwater] => [:percolation], [:x1]),
        SimpleFlux([:infiltration, :percolation, :saturation] => [:outflow]),
        SimpleFlux([:outflow] => [:slowflow, :fastflow]),]

    lfluxes = [
        LagFlux(:slowflow, :x4, LumpedHydro.uh_1_half),
        LagFlux(:fastflow, :x4, LumpedHydro.uh_2_full),
    ]

    dfluxes = [
        StateFlux([:saturation] => [:evap, :percolation], :soilwater, funcs=fluxes)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=vcat(fluxes, lfluxes),
        dfuncs=dfluxes,
        mtk=mtk
    )
end

function FreeWater(; name::Symbol, mtk::Bool=true)

    fluxes = [
        SimpleFlux([:routingstore] => [:recharge], [:x2, :x3, :ω]),
        SimpleFlux([:routingstore] => [:routedflow], [:x3, :γ]),
        SimpleFlux([:routedflow, :recharge, :fastflow_lag] => [:flow])
    ]

    dfluxes = [
        StateFlux([:slowflow_lag] => [:recharge, :routedflow], :routingstore, funcs=fluxes),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=fluxes,
        dfuncs=dfluxes,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true)

    elements = [
        Surface(name=name, mtk=mtk),
        Soil(name=name, mtk=mtk),
        FreeWater(name=name, mtk=mtk)
    ]

    HydroUnit(
        name,
        elements=elements,
    )
end
end