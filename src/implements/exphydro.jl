"""
SoilWaterReservoir in Exp-Hydro
"""
function ExpHydro_SurfElement(; name::Symbol)
    funcs = [
        PetFlux([:temp, :lday]),
        SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DifferFlux([:snowfall], [:melt], :snowwater),
    ]

    HydroElement(
        name=Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function ExpHydro_SoilElement(; name::Symbol)
    funcs = [
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        SurfaceflowFlux([:soilwater], param_names=[:Smax]),
    ]

    dfuncs = [
        DifferFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater)
    ]

    HydroElement(
        name=Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
Inner Route Function in Exphydro
"""
function ExpHydro_RouteElement(; name::Symbol)

    funcs = [
        FlowFlux([:baseflow, :surfaceflow])
    ]

    HydroElement(
        name=name,
        funcs=funcs
    )
end

function ExpHydro_Unit(; name::Symbol)
    elements = [
        ExpHydro_SurfElement(name=name),
        ExpHydro_SoilElement(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function ExpHydro_Node(; name::Symbol)
    HydroNode(
        name,
        units=[ExpHydro_Unit(name=name)],
        routes=namedtuple([name], [ExpHydro_RouteElement(name=name)])
    )
end