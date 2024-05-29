#* baseflow
function expr(eq::HydroEquation{(:soilwater,),(:baseflow,),(:Smax, :Qmax, :f)}; kw...)
    soilwater = first(eq.inputs)
    Smax, Qmax, f = eq.params
    sf = get(kw, :smooth_func, step_func)

    @.[sf(soilwater) * (
        sf(soilwater - Smax) * Qmax +
        sf(soilwater) * sf(Smax - soilwater) * Qmax * exp(-f * (Smax - soilwater))
    )]
end

function expr(eq::HydroEquation{(:routingstore,),(:baseflow,),(:x3, :γ)}; kw...)
    routingstore = first(eq.inputs)
    x3, γ = eq.params

    @.[(x3^(1 - γ)) / (γ - 1) * routingstore^γ]
end

#* totalflow
function expr(eq::HydroEquation{(:baseflow, :surfaceflow),(:totalflow,),()}; kw...)
    baseflow, surfaceflow = eq.inputs

    @.[baseflow + surfaceflow]
end

function expr(eq::HydroEquation{(:routedflow, :recharge, :fastflow),(:totalflow,),()}; kw...)
    routedflow, recharge, fastflow = eq.inputs
    sf = get(kw, :smooth_func, step_func)

    @.[routedflow + sf(fastflow + recharge) * (fastflow + recharge)]
end

function expr(eq::HydroEquation{(:surfaceflow, :baseflow, :interflow),(:totalflow,),()}; kw...)
    surfaceflow, baseflow, interflow = eq.inputs

    @.[surfaceflow + baseflow + interflow]
end

#* surfaceflow
function expr(eq::HydroEquation{(:soilwater,),(:surfaceflow,),(:Smax,)}; kw...)
    soilwater = first(eq.inputs)
    Smax = first(eq.params)
    sf = get(kw, :smooth_func, step_func)

    @.[sf(soilwater) * sf(soilwater - Smax) * (soilwater - Smax)]
end

function expr(eq::HydroEquation{(:surfacerunoff, :prcp),(:surfaceflow,),(:Aim,)}; kw...)
    surfacerunoff, prcp = eq.inputs
    Aim = first(eq.params)

    @.[surfacerunoff + Aim * prcp]
end