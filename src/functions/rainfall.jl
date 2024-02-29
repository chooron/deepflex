function Rainfall(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Rainfall];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=rainfall_func
    )
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:Tmin], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(step_func(input[:Temp] - parameters[:Tmin]) * input[:Prcp])]
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Pet], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(step_func(input[:Prcp] - input[:Pet]) * (input[:Prcp] - input[:Pet]))]
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [input[:Prcp]]
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:tt, :tti], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    tmp_t1 = parameters[:tt] - 0.5 * parameters[:tti]
    tmp_t2 = parameters[:tt] + 0.5 * parameters[:tti]
    [@.(step_func(input[:Temp] - tmp_t2) * input[:Prcp] +
        step_func(tmp_t2 - input[:Temp]) * step_func(input[:Temp] - tmp_t1) * input[:Prcp] * (input[:Temp] - tmp_t1) / parameters[:tti])]
end
