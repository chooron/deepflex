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

function snowfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:tt, :tti], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    tmp_t1 = parameters[:tt] - 0.5 * parameters[:tti]
    tmp_t2 = parameters[:tt] + 0.5 * parameters[:tti]
    [@.(step_func(tmp_t1 - input[:Temp]) * input[:Prcp] +
        step_func(tmp_t2 - input[:Temp]) * step_func(input[:Temp] - tmp_t1) * input[:Prcp] * (tmp_t2 - input[:Temp]) / parameters[:tti])]
end
