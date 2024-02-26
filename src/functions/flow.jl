function Flow(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Flow];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=flow_func
    )
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow, :Surfaceflow], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [input[:Baseflow] .+ input[:Surfaceflow]]
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow, :Recharge, :Q1], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(input[:Baseflow] + step_func(input[:Q1] + input[:Recharge]) * (input[:Q1] + input[:Recharge]))]
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow, :Interflow, :Surfaceflow], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(input[:SurfaceRunoff] + input[:Baseflow] + input[:Interflow])]
end
