function EvapFlux(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:evap;
    param_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=evap_func
    )
end

function evap_func(
    i::namedtuple(:soilwater, :pet),
    p::namedtuple(:Smax),
    sf::Function
)
    @.(sf(i[:soilwater]) * sf(i[:soilwater] - p[:Smax]) * i[:pet] +
    sf(i[:soilwater]) * sf(p[:Smax] - i[:soilwater]) * i[:pet] * (i[:soilwater] / p[:Smax]))
end

function evap_func(
    i::namedtuple(:soilwater, :pet),
    p::namedtuple(:x1),
    sf::Function
)
    @.(i[:pet] * (2 * i[:soilwater] / p[:x1] - (i[:soilwater] / p[:x1])^2))
end

function evap_func(
    i::namedtuple(:soilwater, :pet),
    p::namedtuple(:c, :LM),
    sf::Function
)
    @.(sf(i[:soilwater] - p[:LM]) * i[:pet] +
       sf(p[:LM] - i[:soilwater]) * sf(i[:soilwater] - p[:c] * p[:LM]) * i[:soilwater] / p[:LM] * i[:pet] +
       sf(p[:c] * p[:LM] - i[:soilwater]) * p[:c] * i[:pet])
end

function evap_func(
    i::namedtuple(:soilwater, :pet),
    p::namedtuple(:lp, :fc),
    sf::Function
)
    @.(sf(i[:soilwater] - p[:lp] * p[:fc]) * i[:pet] +
       sf(p[:lp] * p[:fc] - i[:soilwater]) * i[:pet] * i[:soilwater] / (p[:lp] * p[:fc]))
end
