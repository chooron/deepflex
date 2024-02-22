function Surfaceflow(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Surfaceflow],
        parameters,
        surfaceflow_func
    )
end

function surfaceflow_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:Smax], T)
)::gen_namedtuple_type([:Surfaceflow], T) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax = parameters[:Smax]
    (Surfaceflow=@.(step_func(soil_water) * step_func(soil_water - Smax) * (soil_water - Smax)),)
end
