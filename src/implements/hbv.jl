module HBV
"""
SnowWaterReservoir in HyMOD
"""
function SurfElement(; name::Symbol)
    funcs = [
        SnowfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        SimpleFlux([:temp], :refreeze, param_names=[:cfr, :cfmax, :ttm],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:temp]))),
        MeltFlux([:temp], param_names=[:cfmax, :ttm]),
        RainfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        InfiltrationFlux([:snowwater, :liquidwater, :rainfall, :melt], param_names=[:whc]),
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:snowfall, :refreeze], :Out => [:melt]), :snowwater),
        DifferFlux(Dict(:In => [:rainfall, :melt], :Out => [:refreeze, :infiltration]), :liquidwater),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function SoilElement(; name::Symbol)

    funcs = [
        SimpleFlux([:soilwater], :capillary, param_names=[:cflux, :fc],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:soilwater] / p[:fc]))),
        EvapFlux([:soilwater, :pet], param_names=[:lp, :fc]),
        RechargeFlux([:soilwater, :infiltration], param_names=[:fc, :β]),
    ]


    dfuncs = [
        DifferFlux(Dict(:In => [:infiltration, :capillary], :Out => [:evap, :recharge]), :soilwater),
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

function HBVUnit(; name::Symbol)
    elements = [
        SurfaceElement(name=name),
        SoilElement(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function HBVNode(; name::Symbol)
    HydroNode(
        name,
        unit=HBVUnit(name=name),
        route=RouteElement(name=name)
    )
end

end