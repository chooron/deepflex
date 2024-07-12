function expr(eq::HydroEquation{(:soilwater, :infiltration),(:saturation,),(:x1,)}; kw...)
    soilwater, infiltration = eq.inputs
    x1 = first(eq.params)

    [max(0, infiltration * (1 - (soilwater / x1)^2))]
end

function expr(eq::HydroEquation{(:soilwater, :infiltration),(:saturation,),(:Smax, :b)}; kw...)
    soilwater, infiltration = eq.inputs
    Smax, b = eq.params

    [(1 - min(1, max(0, (1 - soilwater / Smax)))^b) * infiltration]
end

function expr(eq::HydroEquation{(:soilwater, :infiltration),(:saturation,),(:Aim, :Wmax, :a, :b)}; kw...)
    soilwater, infiltration = eq.inputs
    Aim, Wmax, a, b = eq.params

    sf = get(kw, :smooth_func, step_func)
    p_i = infiltration .* (1 .- Aim)
    [sf((0.5 - a) - soilwater / Wmax) * (p_i * (abs(0.5 - a)^(1 - b) * abs(soilwater / Wmax)^b)) +
     (sf(soilwater / Wmax - (0.5 - a)) * (p_i * (1 - abs(0.5 + a)^(1 - b) * abs(1 - soilwater / Wmax)^b)))]
end