@reexport module XAJ

using ..LumpedHydro

"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil(; name::Symbol, mtk::Bool=true)
    funcs = [
        SimpleFlux([:prcp], [:infilstration], param_names=[:Aim],
            func=(i, p; kw...) -> @.[(1 - p[:Aim]) * i[:prcp]]),
        SimpleFlux([:infilstration, :soilwater], [:runoff], param_names=[:a, :b, :wmax],
            func=(i, p; kw...) -> @.[i[:infilstration] * (abs(0.5 - p[:a]))^(1 - p[:b]) * (abs(i[:soilwater] / p[:wmax]))^p[:b] +
                                     i[:infilstration] * (1 - abs(0.5 + p[:a]))^(1 - p[:b]) * (abs(1 - i[:soilwater] / p[:wmax]))^p[:b]]),
        EvapFlux([:soilwater, :pet], param_names=[:Smax])
    ]

    dfuncs = [
        StateFlux([:infilstration], [:runoff, :evap], :soilwater, mtk=mtk),
    ]

    HydroElement(
        Symbol(name, :_surface_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

"""
Inner Route Function in Exphydro
"""
function FreeWater(; name::Symbol, mtk=true)

    funcs = [
        FlowFlux([:baseflow, :surfaceflow], :totalflow)
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
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
        SimpleFlux(:totalflow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
    ]

    HydroElement(
        name,
        funcs=funcs,
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