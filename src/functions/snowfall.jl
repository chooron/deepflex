function Snowfall(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Snowfall;
    parameter_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        parameter_names,
        func=snowfall_func
    )
end

function snowfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:Tmin], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(step_func(parameters[:Tmin] - input[:Temp]) * input[:Prcp])
end

function snowfall_func(
    input::gen_namedtuple_type([:Prcp, :Temp], T),
    parameters::gen_namedtuple_type([:tt, :tti], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    tmp_t1 = parameters[:tt] - 0.5 * parameters[:tti]
    tmp_t2 = parameters[:tt] + 0.5 * parameters[:tti]
    @.(step_func(tmp_t1 - input[:Temp]) * input[:Prcp] +
       step_func(tmp_t2 - input[:Temp]) * step_func(input[:Temp] - tmp_t1) * input[:Prcp] * (tmp_t2 - input[:Temp]) / parameters[:tti])
end
