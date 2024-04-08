module GR4J
"""
SnowWaterReservoir in HyMOD
"""
function SurfElement(; name::Symbol)
    funcs = [
        RainfallFlux([:prcp, :pet]),
        SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        InfiltrationFlux([:rainfall])
    ]

    HydroElement(
        name=name,
        funcs=funcs
    )
end

"""
SoilWaterReservoir in HYMOD
"""

"""
SoilWaterReservoir in GR4J
"""
function SoilElement(; name::Symbol)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        EvapFlux([:soilwater, :pet], param_names=[:x1]),
        PercolationFlux([:soilwater], param_names=[:x1]),
        SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        SimpleFlux([:tempflow], :slowflow, param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.9),
        SimpleFlux([:tempflow], :fastflow, param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.1)
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:saturation], :Out => [:evap, :percolation]), :soilwater)
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function RouteElement(; name::Symbol)

    funcs = [
        RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p, sf) -> @.((abs(p[:x3])^(1 - p[:γ])) / (p[:γ] - 1) * (abs(i[:routingstore])^p[:γ]))),
        SimpleFlux([:routedflow, :recharge, :fastflow], :flow,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow]))
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:slowflow, :recharge], :Out => [:routedflow]), :routingstore),
    ]

    lfuncs = [
        LagFlux(:slowflow, :slowflow, lag_func=uh_1_half, param_names=:x4),
        LagFlux(:fastflow, :fastflow, lag_func=uh_2_full, param_names=:x4),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs,
        lfuncs=lfuncs
    )
end

function GR4JUnit(; name::Symbol)
    elements = [
        SurfaceElement(name=name),
        SoilElement(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function GR4JNode(; name::Symbol)
    HydroNode(
        name,
        unit=GR4JUnit(name=name),
        route=RouteElement(name=name)
    )
end

end