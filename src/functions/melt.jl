function Melt(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Melt];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=melt_func
    )
end

function melt_func(
    input::gen_namedtuple_type([:SnowWater, :Temp], T),
    parameters::gen_namedtuple_type([:Tmax, :Df], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    snow_water, temp = input[:SnowWater], input[:Temp]
    Tmax, Df = parameters[:Tmax], parameters[:Df]
    [@.(step_func(temp - Tmax) * step_func(snow_water) * min(snow_water, Df * (temp - Tmax)))]
end

function melt_func(
    input::gen_namedtuple_type([:Temp], T),
    parameters::gen_namedtuple_type([:cfmax, :ttm], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    temp = input[:Temp]
    cfmax, ttm = parameters[:cfmax], parameters[:ttm]
    [@.(step_func(temp - ttm) * (temp - ttm) * cfmax)]
end
