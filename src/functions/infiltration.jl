function InfiltrationFlux(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:infiltration;
    param_names::Vector{Symbol}=Symbol[])

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=infiltration_func
    )
end

function infiltration_func(
    i::namedtuple(:snowwater, :liquidwater, :rainfall, :melt),
    p::namedtuple(:whc),
    sf::Function
)
    @.(sf(i[:liquidwater] - p[:whc] * i[:snowwater]) *
       (i[:rainfall] + i[:melt] + i[:liquidwater] - p[:whc] * i[:snowwater]))
end


function infiltration_func(
    i::namedtuple(:rainfall, :melt),
    p::NamedTuple,
    sf::Function
)
    @.(i[:rainfall] + i[:melt])
end

function infiltration_func(
    i::namedtuple(:rainfall),
    p::NamedTuple,
    sf::Function
)
    i[:rainfall]
end