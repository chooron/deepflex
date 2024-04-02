function SaturationFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:saturation;
    param_names::Vector{Symbol}=Symbol[])

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=saturation_func
    )
end

"""
used for GR4J
"""
function saturation_func(
    i::namedtuple(:soilwater, :infiltration),
    p::namedtuple(:x1),
    sf::Function
)   
    @.(max(0, i[:infiltration] * (1 - (i[:soilwater] / p[:x1])^2)))
end

"""
used in HyMOD
"""
function saturation_func(
    i::namedtuple(:soilwater, :infiltration),
    p::namedtuple(:Smax, :b),
    sf::Function
)
    @.((1 - min(1, max(0, (1 - i[:soilwater] / p[:Smax])))^p[:b]) * i[:infiltration])
end

"""
used in XAJ
"""
function saturation_func(
    i::namedtuple(:soilwater, :infiltration),
    p::namedtuple(:Aim, :Wmax, :a, :b),
    sf::Function
)
    p_i = i[:infiltration] .* (1 .- p[:Aim])
    @.(sf((0.5 - p[:a]) - i[:soilwater] / p[:Wmax]) * (p_i * (abs(0.5 - p[:a])^(1 - p[:b]) * abs(i[:soilwater] / p[:Wmax])^p[:b])) +
       (sf(i[:soilwater] / p[:Wmax] - (0.5 - p[:a])) * (p_i * (1 - (0.5 + p[:a])^(1 - p[:b]) * abs(1 - i[:soilwater] / p[:Wmax])^p[:b]))))
end