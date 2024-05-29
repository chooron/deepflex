include("../../src/LumpedHydro.jl")
using Symbolics
using ModelingToolkit

function test_simpleflux()
    evap_flux = LumpedHydro.SimpleFlux([:soilwater, :pet], [:evap], param_names=[:x1])
    @variables t a(t) b(t) c(t)
    @parameters p1 p2 p3
    custom_flux = LumpedHydro.SimpleFlux([a, b], [c], params=[p1, p2, p3], exprs=[a * p1 + b * p2 + p3])
    custom_flux_func = custom_flux.flux_func
    custom_flux_func((a=1, b=1), (p1=2, p2=3, p3=1))
    # Symbolics.rhss(evap_flux.flux_eqs)
    custom_flux.flux_eqs
end

function test_stateflux()
    #* test state flux
    funcs = [
        LumpedHydro.SimpleFlux([:temp, :lday], [:pet]),
        LumpedHydro.SimpleFlux([:prcp, :temp], [:snowfall], param_names=[:Tmin]),
        LumpedHydro.SimpleFlux([:snowwater, :temp], [:melt], param_names=[:Tmax, :Df]),
        LumpedHydro.SimpleFlux([:prcp, :temp], [:rainfall], param_names=[:Tmin]),
        LumpedHydro.SimpleFlux([:rainfall, :melt], [:infiltration])
    ]

    state_flux = LumpedHydro.StateFlux([:snowfall], [:melt], :snowwater, fluxes=funcs)
    state_flux.state_func((prcp=2.0, temp=3.0, lday=24.0, snowwater=0.0), (Tmin=0.0, Tmax=1.0, Df=3.0))
    state_flux.state_eq
end

function test_neuralflux()
    using Lux
    ann = Lux.Chain(
        Lux.Dense(3 => 16, Lux.tanh),
        Lux.Dense(16 => 1, Lux.leakyrelu)
    )
    nn_flux = LumpedHydro.NeuralFlux([:a, :b, :c], [:d], chain_name=:et, chain=ann)
    tmp_sys = nn_flux.nn_sys
end