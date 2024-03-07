function Pet(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Vector{Symbol}=[:Pet];
    parameters::ComponentVector=ComponentVector(),
    step_func::Function=DEFAULT_SMOOTHER)
    HydroFlux(
        input_names,
        output_names,
        parameters,
        pet_func,
        step_func
    )
end

function pet_func(
    input::gen_namedtuple_type([:Temp, :Lday], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    temp, lday = input[:Temp], input[:Lday]
    @.(29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2))
end