using CSV
using Random
using DataFrames
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OptimizationBBO
using BenchmarkTools
include("../../src/LumpedHydro.jl")

#* parameters in the GR4J model
@parameters x1 = 0.0 [description = "maximum soil moisture storage", unit = "mm"]
@parameters x2 = 0.0 [description = "subsurface water exchange", unit = "mm/d"]
@parameters x3 = 0.0 [description = "routing store depth", unit = "mm"]
@parameters x4 = 0.0 [description = "unit Hydrograph time base", unit = "d"]
#* hydrological flux in the production store
@variables prcp(t) = 0.0 [description = "precipitation", unit = "mm/d"]
@variables ep(t) = 0.0 [description = "potential evaporation", unit = "mm/d"]
@variables soilwater(t) = 0.0 [description = "soil moisture storage", unit = "mm"]
@variables new_soilwater(t) = 0.0 [description = "new soil moisture storage", unit = "mm"]
@variables pn(t) = 0.0 [description = "net precipitation", unit = "mm/d"]
@variables en(t) = 0.0 [description = "net evaporation", unit = "mm/d"]
@variables ps(t) = 0.0 [description = "the fraction of net precipitation redirected to soil moisture", unit = "mm/d"]
@variables es(t) = 0.0 [description = "the fraction of net evaporation subtracted from soil moisture", unit = "mm/d"]
@variables perc(t) = 0.0 [description = "percolation to deeper soil layers", unit = "mm/d"]
@variables pr(t) = 0.0 [description = "combined inflow", unit = "mm/d"]
@variables slowflow(t) = 0.0 [description = "90% ground runoff", unit = "mm/d"]
@variables fastflow(t) = 0.0 [description = "10% direct runoff", unit = "mm/d"]
#* hydrological flux in the unit hydrograph
@variables slowflow_routed(t) = 0.0 [description = "routed ground runoff by uh_1_half", unit = "mm/d"]
@variables fastflow_routed(t) = 0.0 [description = "routed direct runoff by uh_2_full", unit = "mm/d"]
#* hydrological flux in the routing store
@variables routingstore(t) = 0.0 [description = "routing storage", unit = "mm"]
@variables new_routingstore(t) = 0.0 [description = "new routing storage", unit = "mm"]
@variables exch(t) = 0.0 [description = "catchment groundwater exchange", unit = "mm/d"]
@variables routedflow(t) = 0.0 [description = "routed flow in the routing store", unit = "mm/d"]
@variables flow(t) = 0.0 [description = "routed flow in the routing store", unit = "mm/d"]
SimpleFlux = LumpedHydro.SimpleFlux
LagFlux = LumpedHydro.LagFlux
StateFlux = LumpedHydro.StateFlux
HydroElement = LumpedHydro.HydroElement
HydroUnit = LumpedHydro.HydroUnit

#* define the production store
prod_funcs = [
    SimpleFlux([prcp, ep] => [pn, en],
        flux_exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
    SimpleFlux([pn, soilwater] => [ps], [x1],
        flux_exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
    SimpleFlux([en, soilwater] => [es], [x1],
        flux_exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
    SimpleFlux([soilwater] => [perc], [x1],
        flux_exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
    SimpleFlux([pn, ps, perc] => [pr], [x1],
        flux_exprs=[pn - ps + perc]),
    SimpleFlux([pr] => [slowflow, fastflow],
        flux_exprs=[0.9 * pr, 0.1 * pr]),
    SimpleFlux([ps, es, perc, soilwater] => [new_soilwater],
        flux_exprs=[soilwater + ps - es - perc])
]
prod_dfuncs = [StateFlux(soilwater => new_soilwater)]
prod_lfuncs = [
    LagFlux(slowflow => slowflow_routed, x4, LumpedHydro.uh_1_half),
    LagFlux(fastflow => fastflow_routed, x4, LumpedHydro.uh_2_full)
]
prod_ele = HydroElement(:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs, lfuncs=prod_lfuncs)
#* define the routing store
rst_funcs = [
    SimpleFlux([routingstore] => [exch], [x2, x3],
        flux_exprs=[x2 * abs(routingstore / x3)^3.5]),
    SimpleFlux([routingstore, slowflow_routed, exch] => [routedflow], [x3],
        flux_exprs=[x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5]),
    SimpleFlux([routedflow, fastflow_routed, exch] => [flow],
        flux_exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
    SimpleFlux([slowflow_routed, exch, routedflow, routingstore] => [new_routingstore],
        flux_exprs=[routingstore + slowflow_routed + exch - routedflow])
]
rst_dfuncs = [StateFlux(routingstore => new_routingstore)]
rst_ele = HydroElement(:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)
#* define the gr4j model
gr4j_model = HydroUnit(:gr4j, components=[prod_ele, rst_ele])
# load data
df = DataFrame(CSV.File("data/gr4j/sample.csv"));
for col in names(df)[3:end]
    df[ismissing.(df[:, col]), col] .= 0.0
end
prcp_vec = df[!, "prec"]
et_vec = df[!, "pet"]
qobs_vec = df[!, "qobs"]
ts = collect(qobs_vec)

# prepare args
input = (prcp=prcp_vec, ep=et_vec)
pas = ComponentVector(
    params=(x1=320.11, x2=2.42, x3=69.63, x4=1.39),
    initstates=(soilwater=235.97, routingstore=45.47)
)
timeidx = collect(1:length(prcp_vec))
solver = LumpedHydro.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
@benchmark result = gr4j_model(input, pas, timeidx=timeidx, solver=solver); 

# #! set the tunable parameters and constant parameters
# tunable_pas = ComponentVector(params=(x1=320.11, x2=2.42, x3=69.63, x4=1.39))
# const_pas = ComponentVector(initstates=(soilwater=235.97, routingstore=45.47))
# #! set the tunable parameters boundary
# tunable_param_lb = [0.0, -10.0, 1.0, 0.5]
# tunable_param_ub = [2000.0, 10.0, 100.0, 15.0]
# #! prepare flow
# output = (flow=qobs_vec,)
# #! model calibration
# best_pas = LumpedHydro.param_box_optim(
#     gr4j_model,
#     tunable_pas=tunable_pas,
#     const_pas=const_pas,
#     input=input,
#     target=output,
#     timeidx=timeidx,
#     lb=tunable_param_lb,
#     ub=tunable_param_ub,
#     solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
#     maxiters=100,
#     loss_func=LumpedHydro.mse,
#     solver=solver
# )
# #! use the optimized parameters for model simulation
# result_opt = gr4j_model(input, ComponentVector(best_pas; const_pas...), timeidx=timeidx, solver=solver)


# 1 - LumpedHydro.nse(result_opt.flow, df[!, "qobs"])
# 1 - LumpedHydro.nse(result.flow, df[!, "qobs"])

