@reexport module HyMOD
# https://github.com/KMarkert/hymod
using ..LumpedHydro
using ..LumpedHydro.Symbolics: @variables
using ..LumpedHydro.ModelingToolkit: @parameters
using ..LumpedHydro.ModelingToolkit: t_nounits as t
using ..LumpedHydro.ModelingToolkit: Num
"""
SoilWaterReservoir in HYMOD
"""
function SoilStorage(; name::Symbol, mtk::Bool=true)
    @variables prcp = 0.0
    @variables pet = 0.0

    @variables soilwater = 0.0
    @variables tmp_soilwater = 0.0
    @variables new_soilwater = 0.0
    @variables raineff1 = 0.0
    @variables raineff2 = 0.0
    @variables raineff = 0.0
    @variables ct_prev = 0.0
    @variables dummy = 0.0
    @variables evap = 0.0

    @parameters cmax = 0.0
    @parameters bexp = 0.0

    funcs = [
        SimpleFlux([soilwater, prcp, pet] => [ct_prev], [cmax, bexp], flux_exprs=@.[cmax * (1 - abs((1 - ((bexp + 1) * (soilwater) / cmax)))^(1 / (bexp + 1)))]),
        SimpleFlux([prcp, ct_prev] => [raineff1], [cmax], flux_exprs=@.[max((prcp - cmax + ct_prev), 0.0)]),
        SimpleFlux([prcp, ct_prev, raineff1] => [dummy], [cmax], flux_exprs=@.[min(((ct_prev + prcp - raineff1) / cmax), 1.0)]),
        SimpleFlux([dummy] => [tmp_soilwater], [cmax, bexp], flux_exprs=@.[(cmax / (bexp + 1)) * (1 - abs(1 - dummy)^(bexp + 1))]),
        SimpleFlux([soilwater, tmp_soilwater, prcp, raineff1] => [raineff2], Num[], flux_exprs=@.[max(prcp - raineff1 - (tmp_soilwater - soilwater), 0.0)]),
        SimpleFlux([tmp_soilwater, pet] => [evap], [cmax, bexp], flux_exprs=@.[(1 - (((cmax / (bexp + 1)) - tmp_soilwater) / (cmax / (bexp + 1)))) * pet]),
        SimpleFlux([tmp_soilwater, evap] => [new_soilwater], Num[], flux_exprs=@.[max(tmp_soilwater - evap, 0.0)]),
        SimpleFlux([raineff1, raineff2] => [raineff], Num[], flux_exprs=@.[raineff1 + raineff2])
    ]

    dfuncs = [
        StateFlux(new_soilwater => soilwater)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function FreeWaterStorage(; name::Symbol, mtk::Bool=true)
    @variables raineff = 0.0

    @variables slowwater = 0.0
    @variables fastwater1 = 0.0
    @variables fastwater2 = 0.0
    @variables fastwater3 = 0.0

    @variables new_slowwater = 0.0
    @variables new_fastwater1 = 0.0
    @variables new_fastwater2 = 0.0
    @variables new_fastwater3 = 0.0

    @variables slow_q0 = 0.0
    @variables slow_q1 = 0.0
    @variables fast_q0 = 0.0
    @variables fast_q1 = 0.0
    @variables fast_q2 = 0.0
    @variables fast_q3 = 0.0
    @variables flow = 0.0

    @parameters alpha = 0.0
    @parameters ks = 0.0
    @parameters kf = 0.0

    funcs = [
        SimpleFlux([raineff] => [fast_q0, slow_q0], [alpha], flux_exprs=@.[alpha * raineff, (1 - alpha) * raineff]),
        # slow reservoir route
        SimpleFlux([slowwater, slow_q0] => [new_slowwater], [ks], flux_exprs=@.[(1 - ks) * (slowwater + slow_q0)]),
        SimpleFlux([new_slowwater] => [slow_q1], [ks], flux_exprs=@.[(ks / (1 - ks)) * new_slowwater]),
        # fast reservoir route
        SimpleFlux([fastwater1, fast_q0] => [new_fastwater1], [kf], flux_exprs=@.[(1 - kf) * (fastwater1 + fast_q0)]),
        SimpleFlux([new_fastwater1] => [fast_q1], [kf], flux_exprs=@.[(kf / (1 - kf)) * new_fastwater1]),
        SimpleFlux([fastwater2, fast_q1] => [new_fastwater2], [kf], flux_exprs=@.[(1 - kf) * (fastwater2 + fast_q1)]),
        SimpleFlux([new_fastwater2] => [fast_q2], [kf], flux_exprs=@.[(kf / (1 - kf)) * new_fastwater2]),
        SimpleFlux([fastwater3, fast_q2] => [new_fastwater3], [kf], flux_exprs=@.[(1 - kf) * (fastwater3 + fast_q2)]),
        SimpleFlux([new_fastwater3] => [fast_q3], [kf], flux_exprs=@.[(kf / (1 - kf)) * new_fastwater3]),
        # get final output
        SimpleFlux([slow_q1, fast_q3] => [flow], Num[], flux_exprs=@.[max(0.0, slow_q1 + fast_q3)])
    ]

    dfuncs = [
        StateFlux(new_slowwater => slowwater),
        StateFlux(new_fastwater1 => fastwater1),
        StateFlux(new_fastwater2 => fastwater2),
        StateFlux(new_fastwater3 => fastwater3),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function Unit(; name::Symbol, mtk::Bool=true)
    elements = [
        SoilStorage(name=name, mtk=mtk),
        FreeWaterStorage(name=name, mtk=mtk),
    ]

    HydroUnit(
        name,
        components=elements,
    )
end

end