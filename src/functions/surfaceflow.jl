function Surfaceflow(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Vector{Symbol}=[:Surfaceflow];
    parameters::ComponentVector=ComponentVector(),
    step_func::Function=DEFAULT_SMOOTHER)
    HydroFlux(
        input_names,
        output_names,
        parameters,
        surfaceflow_func,
        step_func
    )
end

function surfaceflow_func(
    input::gen_namedtuple_type([:SoilWater], T),
    parameters::gen_namedtuple_type([:Smax], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    soil_water = input[:SoilWater]
    Smax = parameters[:Smax]
    @.(step_func(soil_water) * step_func(soil_water - Smax) * (soil_water - Smax))
end

function surfaceflow_func(
    input::gen_namedtuple_type([:SurfaceRunoff, :Prcp], T),
    parameters::gen_namedtuple_type([:Aim], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    sf, prcp = input[:SurfaceRunoff], input[:Prcp]
    Aim = parameters[:Aim]
    @.(sf + Aim * prcp)
end
