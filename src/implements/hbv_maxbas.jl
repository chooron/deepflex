#* 不考虑融雪径流模块的hbv模型
module HBV_MAXBAS
#* convert from https://github.com/kratzert/RRMPG/blob/master/rrmpg/models/hbvedu_model.py
using ..LumpedHydro
using ..LumpedHydro: step_func
using ..LumpedHydro.Symbolics: @variables
using ..LumpedHydro.ModelingToolkit: @parameters

"""
SoilWaterReservoir in HBV_EDU
"""
function SoilStorage(; name::Symbol)
    @variables soilwater prcpeff prcp pet evap
    @parameters FC Beta PWP

    funcs = [
        SimpleFlux([soilwater, prcp] => [prcpeff], [FC, Beta], flux_exprs=@.[prcp * max(1.0, abs(soilwater / FC))^Beta]),
        SimpleFlux([soilwater, pet] => [evap], [PWP], flux_exprs=@.[(pet * min(1.0, soilwater / PWP))]),
    ]

    dfuncs = [
        StateFlux([prcp] => [prcpeff, evap], soilwater),
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function FreeWaterStorage(; name::Symbol)
    @variables s1 s2 q0 q1 q2 qp prcpeff quz qlz flowuz flowlz
    @parameters L k0 k1 k2 kp conv

    funcs = [
        SimpleFlux([s1] => [q0], [L, k0], flux_exprs=@.[max(0.0, s1 - L) * k0]),
        SimpleFlux([s1] => [q1], [k1], flux_exprs=@.[s1 * k1]),
        SimpleFlux([s1] => [qp], [kp], flux_exprs=@.[s1 * kp]),
        SimpleFlux([s2] => [q2], [k2], flux_exprs=@.[s2 * k2]),
        SimpleFlux([q0] => [quz], flux_exprs=@.[q0]),
        SimpleFlux([q1, q2] => [qlz], flux_exprs=@.[q1 + q2]),
        SimpleFlux([quz] => [flowuz], [conv], flux_exprs=@.[quz * conv]),
        SimpleFlux([qlz] => [flowlz], [conv], flux_exprs=@.[qlz * conv]),
    ]

    dfuncs = [
        StateFlux([prcpeff] => [q0, q1, qp], s1),
        StateFlux([qp] => [q2], s2),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function TriangleRoute(; name::Symbol)
    @variables flowuz_lag flowuz flowlz_lag flowlz
    @parameters lag

    lfuncs = [
        LumpedHydro.LagFlux(flowuz => flowuz_lag, lag, LumpedHydro.uh_3_half),
        # LumpedHydro.LagFlux(flowlz => flowlz_lag, lag, LumpedHydro.uh_3_half),
    ]

    LagElement(
        Symbol(name, :_lag_),
        lfuncs=lfuncs,
    )
end

function OutputElement(; name::Symbol)
    @variables flowlz flowuz_lag flow

    funcs = [
        LumpedHydro.SimpleFlux([flowuz_lag, flowlz] => [flow], flux_exprs=[flowuz_lag + flowlz]),
    ]

    HydroElement(
        Symbol(name, :_lag_),
        funcs=funcs,
    )
end

function Unit(; name::Symbol)
    @info "hello"
    components = [
        SoilStorage(name=name),
        FreeWaterStorage(name=name),
        TriangleRoute(name=name),
        OutputElement(name=name),
    ]

    HydroUnit(
        name,
        components=components
    )
end
end