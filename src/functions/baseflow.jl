function BaseflowFlux(
    i_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:baseflow;
    param_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        i_names,
        output_names,
        param_names=param_names,
        func=baseflow_func,
    )
end

function baseflow_func(
    i::namedtuple(:soilwater),
    p::namedtuple(:Smax, :Qmax, :f),
    sf::Function,
)
    @.(sf(i[:soilwater]) * sf(i[:soilwater] - p[:Smax]) * p[:Qmax] +
    sf(i[:soilwater]) * sf(p[:Smax] - i[:soilwater]) * p[:Qmax] * exp(-p[:f] * (p[:Smax] - i[:soilwater])))
end

function baseflow_func(
    i::namedtuple(:routingstore),
    p::namedtuple(:x3, :γ),
    sf::Function,
)
    @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:routingstore]^p[:γ]))
end

