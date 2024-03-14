function Pet(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Pet;
    param_names::Vector{Symbol}=Symbol[])

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=pet_func
    )
end

function pet_func(
    input::gen_namedtuple_type([:Temp, :Lday], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    temp, lday = input[:Temp], input[:Lday]
    @.(29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2))
end