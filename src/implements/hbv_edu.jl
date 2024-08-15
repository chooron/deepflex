@reexport module HBV_EDU
#* convert from https://github.com/kratzert/RRMPG/blob/master/rrmpg/models/hbvedu_model.py
using ..LumpedHydro
using ..LumpedHydro: step_func
using ..LumpedHydro.Symbolics: @variables
using ..LumpedHydro.ModelingToolkit: @parameters
using ..LumpedHydro.ModelingToolkit: t_nounits as t

"""
SnowWaterReservoir in HBV_EDU
"""
function SnowStorage(; name::Symbol)
    @variables snowwater = 0.0
    @variables liquidwater = 0.0
    @variables prcp = 0.0
    @variables temp = 0.0
    @variables snowfall = 0.0
    @variables melt = 0.0
    @parameters tt = 0.0
    @parameters dd = 0.0

    funcs = [
        SimpleFlux([prcp, temp] => [snowfall], [tt], exprs=@.[step_func(tt - temp) * prcp]),
        SimpleFlux([temp, snowwater] => [melt], [tt, dd], exprs=@.[step_func(temp - tt) * min(snowwater, dd * (temp - tt))]),
        SimpleFlux([prcp, temp, melt] => [liquidwater], [tt], exprs=@.[step_func(temp - tt) * (prcp + melt)]),
    ]

    dfuncs = [
        StateFlux([snowfall] => [melt], snowwater),
    ]

    HydroBucket(
        Symbol(name, :_snow_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

"""
SoilWaterReservoir in HBV_EDU
"""
function SoilStorage(; name::Symbol)
    @variables soilwater = 0.0
    @variables prcpeff = 0.0
    @variables liquidwater = 0.0
    @variables pet = 0.0
    @variables evap = 0.0
    @parameters FC = 0.0
    @parameters Beta = 0.0
    @parameters PWP = 0.0

    funcs = [
        SimpleFlux([soilwater, liquidwater] => [prcpeff], [FC, Beta], exprs=@.[liquidwater * abs(soilwater / FC)^Beta]),
        SimpleFlux([soilwater, pet] => [evap], [PWP], exprs=@.[(pet * min(1.0, soilwater / PWP))]),
    ]

    dfuncs = [
        StateFlux([liquidwater] => [prcpeff, evap], soilwater),
    ]

    HydroBucket(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function FreeWaterStorage(; name::Symbol)
    @variables s1 = 0.0
    @variables s2 = 0.0
    @variables q0 = 0.0
    @variables q1 = 0.0
    @variables q2 = 0.0
    @variables qp = 0.0
    @variables prcpeff = 0.0
    @variables flow = 0.0
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
        StateFlux([prcpeff] => [q0, q1, qp], s1),
        StateFlux([qp] => [q2], s2),
    ]

    HydroBucket(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function Unit(; name::Symbol)
    elements = [
        SnowStorage(name=name),
        SoilStorage(name=name),
        FreeWaterStorage(name=name),
    ]

    HydroUnit(
        name,
        components=elements
    )
end
end