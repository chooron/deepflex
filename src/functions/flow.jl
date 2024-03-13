function Flow(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Flow;
    parameter_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        input_names,
        output_names,
        parameter_names,
        func=flow_func,
    )
end

function flow_func(
    input::gen_namedtuple_type([:BaseFlow, :SurfaceFlow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    input[:BaseFlow] .+ input[:SurfaceFlow]
end

function flow_func(
    input::gen_namedtuple_type([:RoutedFlow, :Recharge, :FastFlow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(input[:RoutedFlow] + step_func(input[:FastFlow] + input[:Recharge]) * (input[:FastFlow] + input[:Recharge]))
end

function flow_func(
    input::gen_namedtuple_type([:BaseFlow, :InterFlow, :SurfaceFlow], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(input[:SurfaceFlow] + input[:BaseFlow] + input[:InterFlow])
end
