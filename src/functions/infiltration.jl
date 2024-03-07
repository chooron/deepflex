function Infiltration(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Vector{Symbol}=[:Infiltration];
    parameters::ComponentVector=ComponentVector(),
    step_func::Function=DEFAULT_SMOOTHER)
    HydroFlux(
        input_names,
        output_names,
        parameters,
        infiltration_func,
        step_func
    )
end

function infiltration_func(
    input::gen_namedtuple_type([:SnowWater, :LiquidWater, :Rainfall, :Melt], T),
    parameters::gen_namedtuple_type([:whc], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(step_func(input[:LiquidWater] - parameters[:whc] * input[:SnowWater]) *
       (input[:Rainfall] + input[:Melt] + input[:LiquidWater] - parameters[:whc] * input[:SnowWater]))
end


function infiltration_func(
    input::gen_namedtuple_type([:Rainfall, :Melt], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(input[:Rainfall] + input[:Melt])
end

function infiltration_func(input::gen_namedtuple_type([:Rainfall], T),
    parameters::NamedTuple,
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    input[:Rainfall]
end