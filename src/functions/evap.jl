function Evap(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:Evap;
    parameters_names::Vector{Symbol}=Symbol[])
    SimpleFlux(
        input_names,
        output_names,
        parameters_names,
        func=evap_func
    )
end

function evap_func(
    input::gen_namedtuple_type([:SoilWater, :Pet], T),
    parameters::gen_namedtuple_type([:Smax], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    Smax = parameters[:Smax]
    @.(step_func(soil_water) * step_func(soil_water - Smax) * pet +
       step_func(soil_water) * step_func(Smax - soil_water) * pet * (soil_water / Smax))
end

function evap_func(
    input::gen_namedtuple_type([:SoilWater, :Pet], T),
    parameters::gen_namedtuple_type([:x1], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    x1 = parameters[:x1]
    @.(pet * (2 * soil_water / x1 - (soil_water / x1)^2))
end

function evap_func(
    input::gen_namedtuple_type([:TensionWater, :Pet], T),
    parameters::gen_namedtuple_type([:c, :LM], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    tension_water, pet = input[:TensionWater], input[:Pet]
    c, LM = parameters[:c], parameters[:LM]
    @.(step_func(tension_water - LM) * pet +
       step_func(LM - tension_water) * step_func(tension_water - c * LM) * tension_water / LM * pet +
       step_func(c * LM - tension_water) * c * pet)
end

function evap_func(
    input::gen_namedtuple_type([:SoilWater, :Pet], T),
    parameters::gen_namedtuple_type([:lp, :fc], T),
    step_func::Function
)::Union{T,Vector{T}} where {T<:Number}
    soil_water, pet = input[:SoilWater], input[:Pet]
    lp, fc = parameters[:lp], parameters[:fc]
    @.(step_func(soil_water - lp * fc) * pet +
       step_func(lp * fc - soil_water) * pet * soil_water / (lp * fc))
end
