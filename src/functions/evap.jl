function expr(eq::HydroEquation{(:soilwater, :pet),(:evap,),(:Smax,)}; kw...)
    soilwater, pet = eq.inputs
    Smax = first(eq.params)

    sf = get(kw, :smooth_func, step_func)

    @.[sf(soilwater) * sf(soilwater - Smax) * pet +
       sf(soilwater) * sf(Smax - soilwater) * pet * (soilwater / Smax)]
end

function expr(eq::HydroEquation{(:soilwater, :pet),(:evap,),(:x1,)}; kw...)
    soilwater, pet = eq.inputs
    x1 = first(eq.params)

    @.[pet * (2 * soilwater / x1 - (soilwater / x1)^2)]
end

function expr(eq::HydroEquation{(:soilwater, :pet),(:evap,),(:c, :LM)}; kw...)
    soilwater, pet = eq.inputs
    c, LM = eq.params

    sf = get(kw, :smooth_func, step_func)
    @.[sf(soilwater - LM) * pet +
       sf(LM - soilwater) * sf(soilwater - c * LM) * soilwater / LM * pet +
       sf(c * LM - soilwater) * c * pet]
end

function expr(eq::HydroEquation{(:soilwater, :pet),(:evap,),(:lp, :fc)}; kw...)
    soilwater, pet = eq.inputs
    lp, fc = eq.params

    sf = get(kw, :smooth_func, step_func)
    @.[sf(soilwater - lp * fc) * pet +
       sf(lp * fc - soilwater) * pet * soilwater / (lp * fc)]
end