function Evap(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Evap],
        parameters,
        evap_func
    )
end

function evap_func(
    input::gen_namedtuple_type([:SoilWater,:Pet], T),
    parameters::gen_namedtuple_type([:Smax], T)
)::gen_namedtuple_type([:Evap], T) where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    Smax = parameters[:Smax]
    (Evap=@.(step_func(soil_water) * step_func(soil_water - Smax) * pet +
             step_func(soil_water) * step_func(Smax - soil_water) * pet * (soil_water / Smax)),)
end

function evap_func(
    input::gen_namedtuple_type([:SoilWater,:Prcp,:Pet], T),
    parameters::gen_namedtuple_type([:x1], T)
)::gen_namedtuple_type([:Evap], T) where {T<:Number}
    soil_water, prcp, pet = input[:SoilWater], input[:Prcp], input[:Pet]
    x1 = parameters[:x1]
    (Evap=@.(step_func(pet - prcp) * (pet - prcp) * (2 * soil_water / x1 - (soil_water / x1)^2)),)
end

