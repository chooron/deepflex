include("../../src/LumpedHydro.jl")

lag_flux = LumpedHydro.LagFlux(:q, :x4, LumpedHydro.uh_2_full)
l = [[1 2 3 2 3 4 1 2]; [1 2 3 2 3 4 1 2]]'
lag_flux.inner_func(l, [1.39])
