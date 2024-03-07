function D_Soilwater(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Vector{Symbol}=[:Soilwater];
    parameters::ComponentVector=ComponentVector(),
    step_func::Function=DEFAULT_SMOOTHER)
    HydroFlux(
        input_names,
        output_names,
        parameters,
        d_soilwater_func,
        step_func
    )
end


function d_soilwater_func(
    input::gen_namedtuple_type([:SoilWater, :Infiltration, :Lday, :Evap, :Flow], T),
    parameters::NamedTuple,
    step_func::Function,
)::Union{T,Vector{T}} where {T<:Number}
    soil_water = input[:SoilWater]
    @.(input[:Infiltration] -
       step_func(soil_water) * input[:Lday] * exp(input[:Evap]) -
       step_func(soil_water) * exp(input[:Flow]))
end

function d_soilwater_func(
    input::gen_namedtuple_type([:SnowWater, :SoilWater, :Rainfall, :Melt, :Lday, :Evap, :Flow], T),
    parameters::NamedTuple,
    step_func::Function,
)::Union{T,Vector{T}} where {T<:Number}
    soil_water = input[:SoilWater]
    @.(relu(sinh(input[:Rainfall])) +
       relu(step_func(input[:SnowWater]) * sinh(input[:Melt])) -
       step_func(soil_water) * input[:Lday] * exp(input[:Evap]) -
       step_func(soil_water) * exp(input[:Flow]))
end