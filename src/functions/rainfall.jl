function expr(eq::HydroEquation{(:prcp, :temp),(:rainfall,),(:Tmin,)}; kw...)
    prcp, temp = eq.inputs
    Tmin = first(eq.params)
    sf = get(kw, :smooth_func, step_func)

    @.[sf(temp - Tmin) * prcp]
end

function expr(eq::HydroEquation{(:prcp, :pet),(:rainfall,),()}; kw...)
    prcp, pet = eq.inputs
    sf = get(kw, :smooth_func, step_func)

    @.[sf(prcp - pet) * (prcp - pet)]
end

function expr(eq::HydroEquation{(:prcp, :pet),(:rainfall,),(:Tmin,)}; kw...)
    prcp, temp = eq.inputs
    tti, tt = eq.params
    sf = get(kw, :smooth_func, step_func)

    tmp_t1 = tt - 0.5 * tti
    tmp_t2 = tt + 0.5 * tti
    @.[sf(temp - tmp_t2) * prcp + sf(tmp_t2 - temp) * sf(temp - tmp_t1) * prcp * (temp - tmp_t1) / tti]
end