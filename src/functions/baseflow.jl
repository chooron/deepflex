function Baseflow(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Baseflow],
        parameters,
        baseflow_func
    )
end

function baseflow_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:Smax, :Qmax, :f], T)
)::gen_namedtuple_type([:Baseflow], T) where {T<:Number}
    soil_water = input[:SoilWater]
    Smax, Qmax, f = parameters[:Smax], parameters[:Qmax], parameters[:f]
    (Baseflow=@.(step_func(soil_water) * step_func(soil_water - Smax) * Qmax +
                 step_func(soil_water) * step_func(Smax - soil_water) * Qmax * exp(-f * (Smax - soil_water))),)
end

