function PetFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:pet;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=pet_func,
        smooth_func=smooth_func
    )
end

function pet_func(
    i::namedtuple(:temp, :lday),
    p::namedtuple();
    kw...
)
    @.(29.8 * i[:lday] * 24 * 0.611 * exp((17.3 * i[:temp]) / (i[:temp] + 237.3)) / (i[:temp] + 273.2))
end

@doc """
hargreaves
"""
function pet_func(
    i::namedtuple(:tmin, :tmax, :datetime),
    p::NamedTuple;
    kw...
)
    tavg = (i[:tmax] .+ i[:tmin]) ./ 2
    b = 2 * pi * (datetime / 365)
    Rav = 1.00011 + 0.034221 * cos(b) + 0.00128 * sin(b) + 0.000719 * cos(2 * b) + 0.000077 * sin(2 * b)
    Ho = ((Gsc * Rav) * 86400) / 1e6
    (0.0023 * Ho * (tmax - tmin)^0.5 * (tavg + 17.8))
end

export PetFlux, pet_func