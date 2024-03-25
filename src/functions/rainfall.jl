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
    i::namedtuple(:prcp, :temp),
    p::namedtuple(:Tmin),
    sf::Function
)
    @.(sf(i[:temp] - p[:Tmin]) * i[:prcp])
end

function rainfall_func(
    i::namedtuple(:prcp, :Pet),
    p::NamedTuple,
    sf::Function
)
    @.(sf(i[:prcp] - i[:pet]) * (i[:prcp] - i[:Pet]))
end

function rainfall_func(
    i::namedtuple(:prcp),
    p::NamedTuple,
    sf::Function
)
    i[:prcp]
end

function rainfall_func(
    i::namedtuple(:prcp, :Temp),
    p::namedtuple(:tt, :tti),
    sf::Function
)
    tmp_t1 = p[:tt] - 0.5 * p[:tti]
    tmp_t2 = p[:tt] + 0.5 * p[:tti]
    @.(sf(i[:Temp] - tmp_t2) * i[:prcp] +
       sf(tmp_t2 - i[:Temp]) * sf(i[:Temp] - tmp_t1) * i[:prcp] * (i[:Temp] - tmp_t1) / p[:tti])
end
