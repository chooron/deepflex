function RainfallFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:rainfall;
    param_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=rainfall_func
    )
end

function rainfall_func(
    i::gen_namedtuple_type([:prcp, :temp], T),
    p::gen_namedtuple_type([:Tmin], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:temp] - p[:Tmin]) * i[:prcp])
end

function rainfall_func(
    i::gen_namedtuple_type([:prcp, :Pet], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:prcp] - i[:pet]) * (i[:prcp] - i[:Pet]))
end

function rainfall_func(
    i::gen_namedtuple_type([:prcp], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    i[:prcp]
end

function rainfall_func(
    i::gen_namedtuple_type([:prcp, :Temp], T),
    p::gen_namedtuple_type([:tt, :tti], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    tmp_t1 = p[:tt] - 0.5 * p[:tti]
    tmp_t2 = p[:tt] + 0.5 * p[:tti]
    @.(sf(i[:Temp] - tmp_t2) * i[:prcp] +
       sf(tmp_t2 - i[:Temp]) * sf(i[:Temp] - tmp_t1) * i[:prcp] * (i[:Temp] - tmp_t1) / p[:tti])
end
