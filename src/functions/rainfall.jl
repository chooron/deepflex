function Rainfall(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Rainfall;
    parameter_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        parameter_names,
        func=rainfall_func
    )
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:Tmin], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(step_func(input[:Temp] - parameters[:Tmin]) * input[:Prcp])
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Pet], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(step_func(input[:Prcp] - input[:Pet]) * (input[:Prcp] - input[:Pet]))
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    input[:Prcp]
end

function rainfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:tt, :tti], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    tmp_t1 = parameters[:tt] - 0.5 * parameters[:tti]
    tmp_t2 = parameters[:tt] + 0.5 * parameters[:tti]
    @.(step_func(input[:Temp] - tmp_t2) * input[:Prcp] +
       step_func(tmp_t2 - input[:Temp]) * step_func(input[:Temp] - tmp_t1) * input[:Prcp] * (input[:Temp] - tmp_t1) / parameters[:tti])
end
