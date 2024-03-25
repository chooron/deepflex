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
    i::namedtuple(:snowwater, :temp),
    p::namedtuple(:Tmax, :Df),
    sf::Function
)
    @.(sf(i[:temp] - p[:Tmax]) * sf(i[:snowwater]) * min(i[:snowwater], p[:Df] * (i[:temp] - p[:Tmax])))
end

function melt_func(
    i::namedtuple(:temp),
    p::namedtuple(:cfmax, :ttm),
    sf::Function
)
    @.(sf(i[:temp] - p[:ttm]) * (i[:temp] - p[:ttm]) * p[:cfmax])
end
