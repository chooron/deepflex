function SurfElement(; name::Symbol,
    funcs::Vector,
    dfuncs::Vector=SimpleFlux[],
    lfuncs::Vector=LagFlux[]
)
    # todo 针对surf element可能有着不同的判断策略

    HydroElement(
        name=Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs,
        lfuncs=lfuncs
    )
end

"""
SnowWaterReservoir in Exp-Hydro
"""
function Surface_ExpHydro(; name::Symbol)
    funcs = [
        PetFlux([:temp, :lday]),
        SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:snowfall], :Out => [:melt]), :snowwater),
    ]

    SurfElement(
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
        RainfallFlux([:prcp, :pet]),
        SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        InfiltrationFlux([:rainfall])
    ]

    SurfElement(
        name=name,
        funcs=funcs
    )
end


function Surface_HBV(; name::Symbol)
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

    SurfElement(
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

    SurfElement(
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

    SurfElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end