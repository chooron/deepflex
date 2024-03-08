function Flow(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Flow;
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=flow_func,
    )
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow, :Surfaceflow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    input[:Baseflow] .+ input[:Surfaceflow]
end

function flow_func(
    input::gen_namedtuple_type([:Routedflow, :Recharge, :Fastflow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(input[:Routedflow] + step_func(input[:Fastflow] + input[:Recharge]) * (input[:Fastflow] + input[:Recharge]))
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow, :Interflow, :Surfaceflow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(input[:SurfaceRunoff] + input[:Baseflow] + input[:Interflow])
end
