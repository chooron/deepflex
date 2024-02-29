function Saturation(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=[:Saturation];
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=saturation_func
    )
end

"""
used for GR4J
"""
function saturation_func(
    input::gen_namedtuple_type([:SoilWater, :Infiltration], T),
    parameters::gen_namedtuple_type([:x1], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(input[:Infiltration] * (1 - (input[:SoilWater] / parameters[:x1])^2))]
end

"""
used in HyMOD
"""
function saturation_func(
    input::gen_namedtuple_type([:SoilWater, :Infiltration], T),
    parameters::gen_namedtuple_type([:Smax, :b], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    soil_water, infiltration = input[:SoilWater], input[:Infiltration]
    Smax, b = parameters[:Smax], parameters[:b]
    @.[(1 - min(1, max(0, (1 - soil_water / Smax)))^b) * infiltration]
end

"""
used in XAJ
"""
function saturation_func(
    input::gen_namedtuple_type([:TensionWater, :Infiltration], T),
    parameters::gen_namedtuple_type([:Aim, :Wmax, :a, :b], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    tension_water, infiltration = input[:TensionWater], input[:Infiltration]
    Aim, Wmax, a, b = parameters[:Aim], parameters[:Wmax], parameters[:a], parameters[:b]
    p_i = infiltration .* (1 .- Aim)
    [@.(step_func((0.5 - a) - tension_water / Wmax) * (p_i * (abs(0.5 - a)^(1 - b) * (tension_water / Wmax)^b)) +
        (step_func(tension_water / Wmax - (0.5 - a)) * (p_i * (1 - (0.5 + a)^(1 - b) * (1 - tension_water / Wmax)^b))))]
end