function expr(eq::HydroEquation{(:prcp, :temp),(:snowfall,),(:Tmin,)}; kw...)
    prcp, temp = eq.inputs
    Tmin = first(eq.params)

    sf = get(kw, :smooth_func, step_func)
    [sf(Tmin - temp) * prcp]
end

function expr(eq::HydroEquation{(:prcp, :temp),(:snowfall,),(:tt, :tti)}; kw...)
    prcp, temp = eq.inputs
    tt, tti = eq.params

    tmp_t1 = tt - 0.5 * tti
    tmp_t2 = tt + 0.5 * tti
    sf = get(kw, :smooth_func, step_func)
    [sf(tmp_t1 - temp) * prcp +
       sf(tmp_t2 - temp) * sf(temp - tmp_t1) * prcp * (tmp_t2 - temp) / tti]
end