function HBV(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}, solver::AbstractSolver) where {T<:Number}
    elements = [
        Surface_HBV(
            name=:sf,
            parameters=parameters[[:tt, :tti, :cfr, :cfmax, :ttm, :whc]],
            init_states=init_states[[:SnowWater, :LiquidWater]]
        ),
        Soil_HBV(
            name=:sl,
            parameters=parameters[[:cflux, :fc, :lp, :k0, :k1, :α, :β, :c]],
            init_states=init_states[[:SoilWater, :UpperZone, :LowerZone]]
        ),
        Lag_HBV(
            name=:lag,
            parameters=parameters[[:maxbas]]
        )
    ]
    build_unit(name=name, elements=elements, solver=solver)
end

function dPL_HBV()
    estimator = [
        
    ]


end