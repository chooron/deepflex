"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_ExpHydro(; name::Symbol)
    funcs = [
        Pet([:temp, :lday]),
        Snowfall([:prcp, :temp], param_names=[:Tmin]),
        Melt([:snowwater, :temp], param_names=[:Tmax, :Df]),
        Rainfall([:prcp, :temp], param_names=[:Tmin]),
        Infiltration([:rainfall, :melt])
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:snowfall], :Out => [:melt]), :snowwater),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_GR4J(; name::Symbol)
    funcs = [
        Rainfall([:prcp, :pet]),
        SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        Infiltration([:rainfall])
    ]

    SimpleElement(
        name=name,
        funcs=funcs
    )
end


function Surface_HBV(; name::Symbol)
    funcs = [
        SnowfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        SimpleFlux([:temp], :refreeze,
            param_names=[:cfr, :cfmax, :ttm],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:temp]))),
        MeltFlux([:temp], param_names=[:cfmax, :ttm]),
        RainfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        InfiltrationFlux([:snowwater, :liquidwater, :rainfall, :melt], param_names=[:whc]),
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:snowfall, :refreeze], :Out => [:melt]), :snowwater),
        DifferFlux(Dict(:In => [:rainfall, :melt], :Out => [:refreeze, :infiltration]), :liquidwater),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function Surface_XAJ(; name::Symbol)
    funcs = [
        SnowfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        SimpleFlux([:temp], :refreeze,
            param_names=[:cfr, :cfmax, :ttm],
            func=(i, p, sf) -> @.(sf(p[:ttm] - i[:temp]) * p[:cfr] * p[:cfmax] * (p[:ttm] - i[:temp]))),
        MeltFlux([:temp], param_names=[:cfmax, :ttm]),
        RainfallFlux([:prcp, :temp], param_names=[:tt, :tti]),
        InfiltrationFlux([:liquidwater, :rainfall, :melt], param_names=[:whc, :sp]),
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:snowfall, :refreeze], :Out => [:melt]), :snowwater),
        DifferFlux(Dict(:In => [:rainfall, :melt], :Out => [:refreeze, :infiltration]), :liquidwater),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


function Surface_M100(; name::Symbol)
    funcs = SimpleFlux[]

    dfuncs = [
        SimpleFlux([:snowwater, :snowfall, :temp, :melt], :snowwater,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(relu(sinh(i[:snowfall]) * sf(i[:temp])) - relu(sf(i[:snowwater]) * sinh(i[:melt])))),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end