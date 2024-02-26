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
    input::gen_namedtuple_type([:SoilWater, :Rainfall], T),
    parameters::gen_namedtuple_type([:x1], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    [@.(input[:Rainfall] * (1 - (input[:SoilWater] / parameters[:x1])^2))]
end

"""
used in HYMOD
"""
function saturation_func(
    input::gen_namedtuple_type([:SoilWater, :Rainfall], T),
    parameters::gen_namedtuple_type([:Smax, :b], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    soil_water, rainfall = input[:Smax], parameters[:Rainfall]
    Smax, b = parameters[:Smax], parameters[:b]
    [@.((1 - min(1, max(0, (1 - soil_water / Smax)))^b) * rainfall)]
end

"""
used in XAJ
"""
function saturation_func(
    input::gen_namedtuple_type([:TensionWater, :Prcp], T),
    parameters::gen_namedtuple_type([:Aim, :Wmax, :a, :b], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    tension_water, prcp = input[:TensionWater], input[:Prcp]
    Aim, Wmax, a, b = parameters[:Aim], parameters[:Wmax], parameters[:a], parameters[:b]
    p_i = prcp .* (1 .- Aim)
    [@.(step_func((0.5 - a) - tension_water / Wmax) * (p_i * (abs(0.5 - a)^(1 - b) * (tension_water / Wmax)^b)) +
        (step_func(tension_water / Wmax - (0.5 - a)) * (p_i * (1 - (0.5 + a)^(1 - b) * (1 - tension_water / Wmax)^b))))]
end

function saturation_func(
    input::gen_namedtuple_type([:FreeWater, :FluxIn], T),
    parameters::gen_namedtuple_type([:Smax, :ex], T)
)::Union{Vector{T},Vector{Vector{T}}} where {T<:Number}
    free_water, flux_in = input[:FreeWater], input[:FluxIn]
    Smax, ex = parameters[:Smax], parameters[:ex]
    tmp_re = @.(step_func(1 - free_water / Smax) * (1 - free_water / Smax))
    [@.(1 - (step_func(1 - tmp_re) * (1 - tmp_re) + step_func(tmp_re - 1) * (tmp_re - 1))^ex * flux_in)]
end