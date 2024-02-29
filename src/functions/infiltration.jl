function Infiltration(
    input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Infiltration];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=infiltration_func
    )
end


function infiltration_func(input::gen_namedtuple_type([:SnowWater, :LiquidWater, :Rainfall, :Melt], T),
    parameters::gen_namedtuple_type([:whc], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    @.[step_func(input[:LiquidWater] - parameters[:whc]* input[:SnowWater]) *
       (input[:Rainfall] + input[:Melt] + input[:LiquidWater] - parameters[:whc] * input[:SnowWater])]
end


function infiltration_func(input::gen_namedtuple_type([:Rainfall, :Melt], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    @.[input[:Rainfall] + input[:Melt]]
end

function infiltration_func(input::gen_namedtuple_type([:Rainfall], T),
    parameters::Nothing=nothing
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [input[:Rainfall]]
end