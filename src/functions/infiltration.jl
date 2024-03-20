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
    i::gen_namedtuple_type([:snowwater, :liquidwater, :rainfall, :melt], T),
    p::gen_namedtuple_type([:whc], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:liquidwater] - p[:whc] * i[:snowwater]) *
       (i[:rainfall] + i[:melt] + i[:liquidwater] - p[:whc] * i[:snowwater]))
end


function infiltration_func(
    i::gen_namedtuple_type([:rainfall, :melt], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(i[:rainfall] + i[:melt])
end

function infiltration_func(
    i::gen_namedtuple_type([:rainfall], T),
    p::NamedTuple,
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    i[:rainfall]
end