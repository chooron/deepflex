function expr(eq::HydroEquation{(:temp, :lday),(:pet,),()}; kw...)
    temp, lday = eq.inputs
    
    @.[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]
end

#* hargreaves
function expr(eq::HydroEquation{(:tmin, :tmax, :datetime),(:pet,),()}; kw...)
    tmin, tmax, datetime = eq.inputs

    tavg = @.((tmax + tmin) / 2)
    b = @.(2 * pi * (datetime / 365))
    Rav = @.(1.00011 + 0.034221 * cos(b) + 0.00128 * sin(b) + 0.000719 * cos(2 * b) + 0.000077 * sin(2 * b))
    Ho = @.(((Gsc * Rav) * 86400) / 1e6)
    @.[0.0023 * Ho * (tmax - tmin)^0.5 * (tavg + 17.8)]
end