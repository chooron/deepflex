using LumpedHydro
using ModelingToolkit

#! construction method 1
perc_flux = SimpleFlux(
    [:soilwater] => [:percolation], [:x1],
    flux_funcs=[(i, p) -> [((p[1])^(-4)) / 4 * ((4 / 9)^(4)) * (i[1]^5)]]
)

#! construction method 2
@variables soilwater = 0.0 [description = "current soil moisture storage"]
@variables percolation = 0.0 [description = "percolation to deeper soil layer"]
@parameters x1 = 0.0 [tunable = true, description = "maximum soil moisture storage"]
perc_flux_v2 = SimpleFlux(
    [soilwater] => [percolation], [x1],
    flux_exprs=[((percolation)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]
)

#! construction method 1
#! common hydrological fluxes (saturation, evaporation, percolation, outflow, etc.)
base_fluxes = [
    SimpleFlux([:soilwater, :infiltration] => [:saturation], [:x1]),
    SimpleFlux([:soilwater, :pet] => [:evap], [:x1]),
    SimpleFlux([:soilwater] => [:percolation], [:x1]),
    SimpleFlux([:infiltration, :percolation, :saturation] => [:outflow]),
    SimpleFlux([:outflow] => [:slowflow, :fastflow])
]
#! state fluxes for soil water within the soil moisture.
dfluxes = [
    StateFlux([:saturation] => [:evap, :percolation], :soilwater, funcs=base_fluxes)
]

#! construction method 2
update_fluxes = [
    SimpleFlux(
        [:soilwater, :saturation, :evap, :percolation] => [:new_soilwater], [:x1],
        flux_funcs=[(i, p) -> [max(p[1], i[1] + i[2] - i[3] - i[4])]]
    )
]
#! state fluxes for soil water within the soil moisture.
dfluxes_v2 = [
    StateFlux([:new_soilwater] => [:soilwater], funcs=vcat(base_fluxes, update_fluxes))
]

#! lag fluxes for slow flow and fast flow.
lfluxes = [
    LagFlux(:slowflow => :slowflow_routed, :x4, LumpedHydro.uh_1_half),
    LagFlux(:fastflow => :fastflow_routed, :x4, LumpedHydro.uh_2_full),
]


#! neural flux for evap and flow.
et_ann = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)
q_ann = Lux.Chain(
    Lux.Dense(2 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)
funcs = [
    StdMeanNormFlux([:snowwater, :soilwater, :prcp, :temp] => [:norm_snw, :norm_slw, :norm_prcp, :norm_temp]),
    NeuralFlux([:norm_snw, :norm_slw, :norm_temp] => [:evap], :etnn => et_ann),
    NeuralFlux([:norm_slw, :norm_prcp] => [:flow], :qnn => q_ann),
]

#! element building
element = HydroElement(:gr4j_soil, funcs=fluxes, dfuncs=dfluxes, lfuncs=lfluxes, mtk=mtk)
#! unit building
unit = HydroUnit(:gr4j, elements=[element])


