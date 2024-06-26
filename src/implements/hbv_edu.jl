@reexport module HBV_EDU
#* convert from https://github.com/kratzert/RRMPG/blob/master/rrmpg/models/hbvedu_model.py
using ..LumpedHydro
using ..LumpedHydro.Symbolics: @variables
using ..LumpedHydro.ModelingToolkit: @parameters
using ..LumpedHydro.ModelingToolkit: t_nounits as t

"""
SnowWaterReservoir in HBV_EDU
"""
function SnowStorage(; name::Symbol, mtk::Bool=true)
    @variables snowwater(t) = 0.0
    @variables liquidwater(t) = 0.0
    @variables prcp(t) = 0.0
    @variables temp(t) = 0.0
    @variables snowfall(t) = 0.0
    @variables melt(t) = 0.0
    @parameters tt = 0.0
    @parameters dd = 0.0

    funcs = [
        SimpleFlux([prcp, temp] => [snowfall], [tt], exprs=@.[step_func(tt - temp) * prcp]),
        SimpleFlux([temp, snowwater] => [melt], [tt, dd], exprs=@.[step_func(temp - tt) * min(snowwater, dd * (temp - tt))]),
        SimpleFlux([prcp, temp, melt] => [liquidwater], [tt], exprs=@.[step_func(temp - tt) * (prcp + melt)]),
    ]

    dfuncs = [
        StateFlux([snowfall] => [melt], snowwater, funcs=funcs),
    ]

    HydroElement(
        Symbol(name, :_snow_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

"""
SoilWaterReservoir in HBV_EDU
"""
function SoilStorage(; name::Symbol, mtk::Bool=true)
    @variables soilwater(t) = 0.0
    @variables prcpeff(t) = 0.0
    @variables liquidwater(t) = 0.0
    @variables pet(t) = 0.0
    @variables evap(t) = 0.0
    @parameters FC = 0.0
    @parameters Beta = 0.0
    @parameters PWP = 0.0

    funcs = [
        SimpleFlux([soilwater, liquidwater] => [prcpeff], [FC, Beta], exprs=@.[liquidwater * abs(soilwater / FC)^Beta]),
        SimpleFlux([soilwater, pet] => [evap], [PWP], exprs=@.[(pet * min(1.0, soilwater / PWP))]),
    ]

    dfuncs = [
        StateFlux([liquidwater] => [prcpeff, evap], soilwater, funcs=funcs),
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function FreeWaterStorage(; name::Symbol, mtk::Bool=true)
    @variables s1(t) = 0.0
    @variables s2(t) = 0.0
    @variables q0(t) = 0.0
    @variables q1(t) = 0.0
    @variables q2(t) = 0.0
    @variables qp(t) = 0.0
    @variables prcpeff(t) = 0.0
    @variables flow(t) = 0.0
    @parameters L = 0.0
    @parameters k0 = 0.0
    @parameters k1 = 0.0
    @parameters k2 = 0.0
    @parameters kp = 0.0
    @parameters PWP = 0.0

    funcs = [
        SimpleFlux([s1] => [q0], [L, k0], exprs=@.[max(0.0, s1 - L) * k0]),
        SimpleFlux([s1] => [q1], [k1], exprs=@.[s1 * k1]),
        SimpleFlux([s1] => [qp], [kp], exprs=@.[s1 * kp]),
        SimpleFlux([s2] => [q2], [k2], exprs=@.[s2 * k2]),
        SimpleFlux([q0, q1, q2] => [flow], exprs=@.[q0 + q1 + q2]),
    ]

    dfuncs = [
        StateFlux([prcpeff] => [q0, q1, qp], s1, funcs=funcs),
        StateFlux([qp] => [q2], s2, funcs=funcs),
    ]

    HydroElement(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true, )
    elements = [
        SnowStorage(name=name, mtk=mtk),
        SoilStorage(name=name, mtk=mtk),
        FreeWaterStorage(name=name, mtk=mtk),
    ]

    HydroUnit(
        name,
        elements=elements
    )
end
end