function BaseFlow(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:BaseFlow;
    param_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=baseflow_func,
    )
end

function baseflow_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:Smax, :Qmax, :f], T),
    step_func::Function,
)::Union{T,Vector{T}} where {T<:Number}
    soil_water = input[:SoilWater]
    Smax, Qmax, f = parameters[:Smax], parameters[:Qmax], parameters[:f]
    @.(step_func(soil_water) * step_func(soil_water - Smax) * Qmax +
       step_func(soil_water) * step_func(Smax - soil_water) * Qmax * exp(-f * (Smax - soil_water)))
end

function baseflow_func(
    input::gen_namedtuple_type([:RoutingStore], T),
    parameters::gen_namedtuple_type([:x3, :γ], T),
    step_func::Function,
)::Union{T,Vector{T}} where {T<:Number}
    routing_store = input[:RoutingStore]
    x3, γ = parameters[:x3], parameters[:γ]
    @.((x3^(1 - γ)) / (γ - 1) * (routing_store^γ))
end

