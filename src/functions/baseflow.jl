function Baseflow(
    input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Baseflow];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=baseflow_func
    )
end

function baseflow_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:Smax, :Qmax, :f], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    soil_water = input[:SoilWater]
    Smax, Qmax, f = parameters[:Smax], parameters[:Qmax], parameters[:f]
    [@.(step_func(soil_water) * step_func(soil_water - Smax) * Qmax + step_func(soil_water) * step_func(Smax - soil_water) * Qmax * exp(-f * (Smax - soil_water)))]
end

function baseflow_func(
    input::gen_namedtuple_type([:RoutingStore], T),
    parameters::gen_namedtuple_type([:x3, :γ], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    routing_store = input[:RoutingStore]
    x3, γ = parameters[:x3], parameters[:γ]
    [@.((x3^(1 - γ)) / (γ - 1) * (routing_store^γ))]
end

