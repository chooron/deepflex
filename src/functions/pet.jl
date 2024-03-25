function PetFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:pet;
    param_names::Vector{Symbol}=Symbol[])

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=pet_func
    )
end

function pet_func(
    i::namedtuple(:temp, :lday),
    p::namedtuple(),
    sf::Function
)
    @.(29.8 * i[:lday] * 24 * 0.611 * exp((17.3 * i[:temp]) / (i[:temp] + 237.3)) / (i[:temp] + 273.2))
end