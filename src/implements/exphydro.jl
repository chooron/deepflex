@reexport module ExpHydro

using ..DeepFlex

"""
SoilWaterReservoir in Exp-Hydro
"""
function SurfElement(; name::Symbol)
    funcs = [
        DeepFlex.PetFlux([:temp, :lday]),
        DeepFlex.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        DeepFlex.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:snowfall], [:melt], :snowwater),
    ]

    DeepFlex.HydroElement(
        name=Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs
    )
end

# """
# SoilWaterReservoir in Exp-Hydro
# """
# function SoilElement(; name::Symbol)
#     funcs = [
#         EvapFlux([:soilwater, :pet], param_names=[:Smax]),
#         BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
#         SurfaceflowFlux([:soilwater], param_names=[:Smax]),
#     ]

#     dfuncs = [
#         DifferFlux([:infiltration], [:evap, :baseflow, :surfaceflow], :soilwater)
#     ]

#     HydroElement(
#         name=Symbol(name, :_soil_),
#         funcs=funcs,
#         dfuncs=dfuncs
#     )
# end

# """
# Inner Route Function in Exphydro
# """
# function ZoneElement(; name::Symbol)

#     funcs = [
#         FlowFlux([:baseflow, :surfaceflow])
#     ]

#     HydroElement(
#         name=name,
#         funcs=funcs
#     )
# end

# function Route(; name::Symbol)

#     funcs = [
#         SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
#     ]

#     RouteElement(
#         name=name,
#         funcs=funcs
#     )
# end

# function Unit(; name::Symbol)
#     elements = [
#         SurfElement(name=name),
#         SoilElement(name=name),
#         ZoneElement(name=name),
#     ]
#     HydroUnit(name, elements=elements)
# end

# function Node(; name::Symbol)
#     HydroNode(
#         name,
#         units=[Unit(name=name)],
#         routes=namedtuple([name], [Route(name=name)])
#     )
# end
end