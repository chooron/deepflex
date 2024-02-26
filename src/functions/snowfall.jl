function Snowfall(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Snowfall];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=snowfall_func
    )
end

function snowfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:Tmin], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(step_func(parameters[:Tmin] - input[:Temp]) * input[:Prcp])]
end
