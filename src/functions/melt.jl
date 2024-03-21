function MeltFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:melt;
    param_names::Vector{Symbol}=Symbol[])

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=melt_func
    )
end

function melt_func(
    i::gen_namedtuple_type([:snowwater, :temp], T),
    p::gen_namedtuple_type([:Tmax, :Df], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:temp] - p[:Tmax]) * sf(i[:snowwater]) * min(i[:snowwater], p[:Df] * (i[:temp] - p[:Tmax])))
end

function melt_func(
    i::gen_namedtuple_type([:temp], T),
    p::gen_namedtuple_type([:cfmax, :ttm], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:temp] - p[:ttm]) * (i[:temp] - p[:ttm]) * p[:cfmax])
end
