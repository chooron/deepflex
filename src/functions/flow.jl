#* baseflow
function expr(eq::HydroEquation{(:soilwater,),(:baseflow,),(:Smax, :Qmax, :f)}; kw...)
    soilwater = first(eq.inputs)
    Smax, Qmax, f = eq.params
    sf = get(kw, :smooth_func, step_func)

    @.[sf(soilwater) * (
        sf(soilwater - Smax) * Qmax +
        sf(Smax - soilwater) * Qmax * exp(-f * (Smax - soilwater))
    )]
end

function expr(eq::HydroEquation{(:routingstore,),(:baseflow,),(:x3, :γ)}; kw...)
    routingstore = first(eq.inputs)
    x3, γ = eq.params

    @.[(x3^(1 - γ)) / (γ - 1) * routingstore^γ]
end

#* flow
function expr(eq::HydroEquation{(:baseflow, :surfaceflow),(:flow,),()}; kw...)
    baseflow, surfaceflow = eq.inputs

    @.[baseflow + surfaceflow]
end

function expr(eq::HydroEquation{(:routedflow, :recharge, :fastflow),(:flow,),()}; kw...)
    routedflow, recharge, fastflow = eq.inputs
    sf = get(kw, :smooth_func, step_func)

    @.[routedflow + sf(fastflow + recharge) * (fastflow + recharge)]
end

function expr(eq::HydroEquation{(:surfaceflow, :baseflow, :interflow),(:flow,),()}; kw...)
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

# for gr4j
function expr(eq::HydroEquation{(:infiltration, :percolation, :saturation),(:outflow,),()}; kw...)
    infiltration, percolation, saturation = eq.inputs
    @.[infiltration - saturation + percolation]
end

function expr(eq::HydroEquation{(:outflow,),(:slowflow, :fastflow),()}; kw...)
    outflow = first(eq.inputs)
    @.[outflow * 0.9, outflow * 0.1]
end

function expr(eq::HydroEquation{(:routingstore, :recharge, :slowflow_routed),(:routedflow,),(:x3, :γ)}; kw...)
    routingstore, recharge, slowflow_lag = eq.inputs
    rgt = @.(routingstore + slowflow_lag + recharge)
    x3, γ = eq.params
    @.[((abs(x3)^(1 - γ)) / (γ - 1) * (abs(rgt)^γ))]
    # @.[rgt * (1.0 - (1.0 + (rgt / x3)^4)^(-1 / 4))]
end

function expr(eq::HydroEquation{(:routedflow, :recharge, :fastflow_routed),(:flow,),()}; kw...)
    routedflow, recharge, fastflow_lag = eq.inputs
    @.[routedflow + max(0.0, recharge + fastflow_lag)]
end