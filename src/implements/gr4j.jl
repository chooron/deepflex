@reexport module GR4J

using ..HydroModels

"""
SoilWaterReservoir in GR4J
"""
function Surface(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:prcp, :pet] => [:rainfall]),
        SimpleFlux([:rainfall] => [:infiltration])
    ]

    HydroBucket(
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
        LagFlux(:slowflow => :slowflow_routed, :x4, HydroEquations.uh_1_half),
        LagFlux(:fastflow => :fastflow_routed, :x4, HydroEquations.uh_2_full),
    ]

    dfluxes = [
        StateFlux([:saturation] => [:evap, :percolation], :soilwater)
    ]

    HydroBucket(
        name=Symbol(name, :_soil_),
        funcs=fluxes,
        dfuncs=dfluxes,
        lfuncs=lfluxes,
    )
end

function FreeWater(; name::Symbol)

    fluxes = [
        SimpleFlux([:routingstore] => [:recharge], [:x2, :x3, :ω]),
        SimpleFlux([:routingstore, :recharge, :slowflow_routed] => [:routedflow], [:x3, :γ]),
        SimpleFlux([:routedflow, :recharge, :fastflow_routed] => [:flow])
    ]

    dfluxes = [
        StateFlux([:slowflow_routed, :recharge] => [:routedflow], :routingstore),
    ]

    HydroBucket(
        name=Symbol(name, :_zone_),
        funcs=fluxes,
        dfuncs=dfluxes,
    )
end

function Model(; name::Symbol)

    elements = [
        Surface(name=name),
        Soil(name=name),
        FreeWater(name=name)
    ]

    HydroModel(
        name,
        components=elements,
    )
end
end