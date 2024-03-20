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
    i::gen_namedtuple_type([:soilwater], T),
    p::gen_namedtuple_type([:Smax, :Qmax, :f], T),
    sf::Function,
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:soilwater]) * sf(i[:soilwater] - p[:Smax]) * p[:Qmax] +
    sf(i[:soilwater]) * sf(p[:Smax] - i[:soilwater]) * p[:Qmax] * exp(-p[:f] * (p[:Smax] - i[:soilwater])))
end

function baseflow_func(
    i::gen_namedtuple_type([:routingstore], T),
    p::gen_namedtuple_type([:x3, :γ], T),
    sf::Function,
)::Union{T,Vector{T}} where {T<:Number}
    @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:routingstore]^p[:γ]))
end

