"""
Exp-Hydro model
"""
function ExpHydro(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}, solver::AbstractSolver=nothing) where {T<:Number}
    elements = [
        Surface_ExpHydro(
            name=:sf,
            parameters=parameters[[:Tmin, :Tmax, :Df]],
            init_states=init_states[[:SnowWater]]
        ),
        Soil_ExpHydro(
            name=:sl,
            parameters=parameters[[:Smax, :Qmax, :f]],
            init_states=init_states[[:SoilWater]]
        ),
        Routing_ExpHydro(name=:rt)
    ]
    build_unit(name=name, elements=elements, solver=solver)
end
