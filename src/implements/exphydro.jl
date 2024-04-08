module ExpHydro
"""
SoilWaterReservoir in Exp-Hydro
"""
function SurfaceElement(; name::Symbol)
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

    HydroElement(
        name=Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SoilWaterReservoir in Exp-Hydro
"""
function SoilElement(; name::Symbol)
    funcs = [
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        SurfaceflowFlux([:soilwater], param_names=[:Smax]),
        FlowFlux([:baseflow, :surfaceflow])
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:infiltration], :Out => [:evap, :flow]), :soilwater)
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
function RouteElement(; name::Symbol)

    funcs = [
        FlowFlux([:baseflow, :surfaceflow])
    ]

    HydroElement(
        name=name,
        funcs=funcs
    )
end

function ExphydroUnit(; name::Symbol)
    elements = [
        SurfaceElement(name=name),
        SoilElement(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function ExphydroNode(; name::Symbol)
    HydroNode(
        name,
        unit=ExphydroUnit(name=name),
        route=RouteElement(name=name)
    )
end


end
