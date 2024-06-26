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
    @variables prcp(t) = 0.0
    @variables pet(t) = 0.0

    @variables soilwater(t) = 0.0
    @variables tmp_soilwater(t) = 0.0
    @variables new_soilwater(t) = 0.0
    @variables raineff1(t) = 0.0
    @variables raineff2(t) = 0.0
    @variables raineff(t) = 0.0
    @variables ct_prev(t) = 0.0
    @variables dummy(t) = 0.0
    @variables evap(t) = 0.0

    @parameters cmax = 0.0
    @parameters bexp = 0.0

    funcs = [
        SimpleFlux([soilwater, prcp, pet] => [ct_prev], [cmax, bexp], exprs=@.[cmax * (1 - abs((1 - ((bexp + 1) * (soilwater) / cmax)))^(1 / (bexp + 1)))]),
        SimpleFlux([prcp, ct_prev] => [raineff1], [cmax], exprs=@.[max((prcp - cmax + ct_prev), 0.0)]),
        SimpleFlux([prcp, ct_prev, raineff1] => [dummy], [cmax], exprs=@.[min(((ct_prev + prcp - raineff1) / cmax), 1.0)]),
        SimpleFlux([dummy] => [tmp_soilwater], [cmax, bexp], exprs=@.[(cmax / (bexp + 1)) * (1 - abs(1 - dummy)^(bexp + 1))]),
        SimpleFlux([soilwater, tmp_soilwater, prcp, raineff1] => [raineff2], Num[], exprs=@.[max(prcp - raineff1 - (tmp_soilwater - soilwater), 0.0)]),
        SimpleFlux([tmp_soilwater, pet] => [evap], [cmax, bexp], exprs=@.[(1 - (((cmax / (bexp + 1)) - tmp_soilwater) / (cmax / (bexp + 1)))) * pet]),
        SimpleFlux([tmp_soilwater, evap] => [new_soilwater], Num[], exprs=@.[max(tmp_soilwater - evap, 0.0)]),
        SimpleFlux([raineff1, raineff2] => [raineff], Num[], exprs=@.[raineff1 + raineff2])
    ]

    dfuncs = [
        StateFlux(new_soilwater => soilwater, funcs=funcs)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk,
    )
end

function FreeWaterStorage(; name::Symbol, mtk::Bool=true)
    @variables raineff(t) = 0.0

    @variables slowwater(t) = 0.0
    @variables fastwater1(t) = 0.0
    @variables fastwater2(t) = 0.0
    @variables fastwater3(t) = 0.0

    @variables new_slowwater(t) = 0.0
    @variables new_fastwater1(t) = 0.0
    @variables new_fastwater2(t) = 0.0
    @variables new_fastwater3(t) = 0.0

    @variables slow_q0(t) = 0.0
    @variables slow_q1(t) = 0.0
    @variables fast_q0(t) = 0.0
    @variables fast_q1(t) = 0.0
    @variables fast_q2(t) = 0.0
    @variables fast_q3(t) = 0.0
    @variables flow(t) = 0.0

    @parameters alpha = 0.0
    @parameters ks = 0.0
    @parameters kf = 0.0

    funcs = [
        SimpleFlux([raineff] => [fast_q0, slow_q0], [alpha], exprs=@.[alpha * raineff, (1 - alpha) * raineff]),
        # slow reservoir route
        SimpleFlux([slowwater, slow_q0] => [new_slowwater], [ks], exprs=@.[(1 - ks) * (slowwater + slow_q0)]),
        SimpleFlux([new_slowwater] => [slow_q1], [ks], exprs=@.[(ks / (1 - ks)) * new_slowwater]),
        # fast reservoir route
        SimpleFlux([fastwater1, fast_q0] => [new_fastwater1], [kf], exprs=@.[(1 - kf) * (fastwater1 + fast_q0)]),
        SimpleFlux([new_fastwater1] => [fast_q1], [kf], exprs=@.[(kf / (1 - kf)) * new_fastwater1]),
        SimpleFlux([fastwater2, fast_q1] => [new_fastwater2], [kf], exprs=@.[(1 - kf) * (fastwater2 + fast_q1)]),
        SimpleFlux([new_fastwater2] => [fast_q2], [kf], exprs=@.[(kf / (1 - kf)) * new_fastwater2]),
        SimpleFlux([fastwater3, fast_q2] => [new_fastwater3], [kf], exprs=@.[(1 - kf) * (fastwater3 + fast_q2)]),
        SimpleFlux([new_fastwater3] => [fast_q3], [kf], exprs=@.[(kf / (1 - kf)) * new_fastwater3]),
        # get final output
        SimpleFlux([slow_q1, fast_q3] => [flow], Num[], exprs=@.[max(0.0, slow_q1 + fast_q3)])
    ]

    dfuncs = [
        StateFlux(new_slowwater => slowwater, funcs=funcs),
        StateFlux(new_fastwater1 => fastwater1, funcs=funcs),
        StateFlux(new_fastwater2 => fastwater2, funcs=funcs),
        StateFlux(new_fastwater3 => fastwater3, funcs=funcs),
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
        elements=elements,
    )
end

end