function expr(eq::HydroEquation{(:routingstore,),(:recharge,),(:x2, :x3, :ω)}; kw...)
    routingstore = first(eq.inputs)
    x2, x3, ω = eq.params

    [x2 / (abs(x3)^ω) * abs(routingstore)^ω]
end

function expr(eq::HydroEquation{(:soilwater, :infiltration),(:recharge,),(:fc, :β)}; kw...)
    soilwater, infiltration = eq.inputs
    fc, β = eq.params

    [infiltration * (soilwater / fc)^β]
end