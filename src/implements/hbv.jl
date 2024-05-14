@reexport module HBV

using ..LumpedHydro
"""
SnowWaterReservoir in HyMOD
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        SnowfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        SimpleFlux([:temp], :refreeze, param_names=[:cfr, :cfmax, :ttm],
            func=(i, p; kw...) -> begin
                sf = kw[:smooth_func]
                @.(sf(p[:ttm] - i[:temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:temp]))
            end),
        MeltFlux([:temp], param_names=[:cfmax, :ttm]),
        RainfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        InfiltrationFlux([:snowwater, :liquidwater, :rainfall, :melt], param_names=[:whc]),
    ]

    dfuncs = [
        StateFlux([:snowfall, :refreeze], [:melt], :snowwater),
        StateFlux([:rainfall, :melt], [:refreeze, :infiltration], :liquidwater),
    ]

    HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Soil(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:soilwater], :capillary, param_names=[:cflux, :fc],
            func=(i, p; kw...) -> @.(p[:cflux] * (1 - i[:soilwater] / p[:fc]))),
        EvapFlux([:soilwater, :pet], param_names=[:lp, :fc]),
        RechargeFlux([:soilwater, :infiltration], param_names=[:fc, :β]),
    ]

    dfuncs = [
        StateFlux([:infiltration, :capillary], [:evap, :recharge], :soilwater),
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function FreeWater(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux([:upperzone], [:q0], param_names=[:k0, :α],
            func=(i, p; kw...) -> @.(p[:k0] .* (i[:upperzone]^(1 + p[:α])))),
        SimpleFlux([:lowerzone], [:q1], param_names=[:k1],
            func=(i, p; kw...) -> @.(p[:k1] .* i[:lowerzone])),
        SimpleFlux(Symbol[], :percolation, param_names=[:c],
            func=(i, p; kw...) -> p[:c]),
        SimpleFlux([:q0, :q1], :totalflow, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:q0] + i[:q1]))
    ]

    dfuncs = [
        StateFlux([:recharge], [:capillary, :q0, :percolation], :upperzone),
        StateFlux([:percolation], [:q1], :lowerzone),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroUnit(
        name,
        elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk), FreeWater(name=name, mtk=mtk)],
        step=step,
    )
end

function Route(; name::Symbol)

    funcs = [
        LagFlux(:totalflow, :flow, lag_func=LumpedHydro.uh_4_full, lag_time=:maxbas),
    ]

    HydroElement(
        name,
        funcs=funcs,
        mtk=false
    )
end

function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroNode(
        name,
        units=[Unit(name=name, mtk=mtk, step=step)],
        routes=[Route(name=name)],
    )
end

end