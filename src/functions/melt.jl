function expr(eq::HydroEquation{(:snowwater, :temp),(:melt,),(:Tmax, :Df)};kw...)
    snowwater, temp = eq.inputs
    Tmax, Df = eq.params
    sf = get(kw, :smooth_func, step_func)

    [sf(temp - Tmax) * sf(snowwater) * min(snowwater, Df * (temp - Tmax))]
end

function expr(eq::HydroEquation{(:temp,),(:melt,),(:cfmax, :ttm)};kw...)
    temp = first(eq.inputs)
    cfmax, ttm = eq.params

    sf = get(kw, :smooth_func, step_func)
    [sf(temp - ttm) * (temp - ttm) * cfmax]
end