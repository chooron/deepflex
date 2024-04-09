"""
SnowWaterReservoir in HyMOD
"""
function HyMOD_SurfElement(; name::Symbol)
    funcs = [
        RainfallFlux([:prcp, :pet]),
        SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        InfiltrationFlux([:rainfall])
    ]

    HydroElement(
        name=name,
        funcs=funcs
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function HyMOD_SoilElement(; name::Symbol)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:Smax, :b]),
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        SimpleFlux([:saturation], :fastflow, param_names=[:a],
            func=(i, p, sf) -> @.(i[:saturation] * (1 - p[:a]))),
        SimpleFlux([:saturation], :slowflow, param_names=[:a],
            func=(i, p, sf) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:saturation], :Out => [:evap, :saturation]), :soilwater)
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


function HyMOD_RouteElement(; name::Symbol)

    funcs = [
        SimpleFlux([:fr1], :qf1, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr1]),
        SimpleFlux([:fr2], :qf2, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr2]),
        SimpleFlux([:fr3], :qf3, param_names=[:kf], func=(i, p, sf) -> p[:kf] .* i[:fr3]),
        SimpleFlux([:sr], :qs, param_names=[:ks], func=(i, p, sf) -> p[:ks] .* i[:sr]),
        SimpleFlux([:qs, :qf3], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:qs] .+ i[:qf3]),
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:fastflow], :Out => [:qf1]), :fr1),
        DifferFlux(Dict(:In => [:qf1], :Out => [:qf2]), :fr2),
        DifferFlux(Dict(:In => [:qf2], :Out => [:qf3]), :fr3),
        DifferFlux(Dict(:In => [:slowflow], :Out => [:qs]), :sr),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function HyMOD_Unit(; name::Symbol)
    elements = [
        HyMOD_SurfElement(name=name),
        HyMOD_SoilElement(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function HyMOD_Node(; name::Symbol)
    HydroNode(
        name,
        unit=HyMOD_Unit(name=name),
        route=HyMOD_RouteElement(name=name)
    )
end