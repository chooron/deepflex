using CSV
using Random
using DataFrames
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using OptimizationBBO
using BenchmarkTools
using HydroErrors
include("../../src/HydroModels.jl")

#* parameters in the GR4J model
@parameters x1 = 0.0 [description = "maximum soil moisture storage", unit = "mm"]
@parameters x2 = 0.0 [description = "subsurface water exchange", unit = "mm/d"]
@parameters x3 = 0.0 [description = "routing store depth", unit = "mm"]
@parameters x4 = 0.0 [description = "unit Hydrograph time base", unit = "d"]
#* hydrological flux in the production store
@variables prcp = 0.0 [description = "precipitation", unit = "mm/d"]
@variables ep = 0.0 [description = "potential evaporation", unit = "mm/d"]
@variables soilwater = 0.0 [description = "soil moisture storage", unit = "mm"]
@variables new_soilwater = 0.0 [description = "new soil moisture storage", unit = "mm"]
@variables pn = 0.0 [description = "net precipitation", unit = "mm/d"]
@variables en = 0.0 [description = "net evaporation", unit = "mm/d"]
@variables ps = 0.0 [description = "the fraction of net precipitation redirected to soil moisture", unit = "mm/d"]
@variables es = 0.0 [description = "the fraction of net evaporation subtracted from soil moisture", unit = "mm/d"]
@variables perc = 0.0 [description = "percolation to deeper soil layers", unit = "mm/d"]
@variables pr = 0.0 [description = "combined inflow", unit = "mm/d"]
@variables slowflow = 0.0 [description = "90% ground runoff", unit = "mm/d"]
@variables fastflow = 0.0 [description = "10% direct runoff", unit = "mm/d"]
#* hydrological flux in the unit hydrograph
@variables slowflow_routed = 0.0 [description = "routed ground runoff by uh_1_half", unit = "mm/d"]
@variables fastflow_routed = 0.0 [description = "routed direct runoff by uh_2_full", unit = "mm/d"]
#* hydrological flux in the routing store
@variables routingstore = 0.0 [description = "routing storage", unit = "mm"]
@variables new_routingstore = 0.0 [description = "new routing storage", unit = "mm"]
@variables exch = 0.0 [description = "catchment groundwater exchange", unit = "mm/d"]
@variables routedflow = 0.0 [description = "routed flow in the routing store", unit = "mm/d"]
@variables flow = 0.0 [description = "routed flow in the routing store", unit = "mm/d"]
SimpleFlux = HydroModels.SimpleFlux
UnitHydroRouteFlux = HydroModels.UnitHydroRouteFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel

#* define the production store
prod_funcs = [
    SimpleFlux([prcp, ep] => [pn, en],
        exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
    SimpleFlux([pn, soilwater] => [ps], [x1],
        exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
    SimpleFlux([en, soilwater] => [es], [x1],
        exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
    SimpleFlux([soilwater] => [perc], [x1],
        exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
    SimpleFlux([pn, ps, perc] => [pr], [x1],
        exprs=[pn - ps + perc]),
    SimpleFlux([pr] => [slowflow, fastflow],
        exprs=[0.9 * pr, 0.1 * pr]),
    SimpleFlux([ps, es, perc, soilwater] => [new_soilwater],
        exprs=[soilwater + ps - es - perc])
]
prod_dfuncs = [StateFlux(soilwater => new_soilwater)]

uh_flux_1 = UnitHydroRouteFlux(slowflow, x4, HydroModels.uh_1_half)
uh_flux_2 = UnitHydroRouteFlux(fastflow, x4, HydroModels.uh_2_full)

prod_ele = HydroBucket(:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs)
#* define the routing store
rst_funcs = [
    SimpleFlux([routingstore] => [exch], [x2, x3],
        exprs=[x2 * abs(routingstore / x3)^3.5]),
    SimpleFlux([routingstore, slowflow_routed, exch] => [routedflow], [x3],
        exprs=[x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5]),
    SimpleFlux([routedflow, fastflow_routed, exch] => [flow],
        exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
    SimpleFlux([slowflow_routed, exch, routedflow, routingstore] => [new_routingstore],
        exprs=[routingstore + slowflow_routed + exch - routedflow])
]
rst_dfuncs = [StateFlux(routingstore => new_routingstore)]
rst_ele = HydroBucket(:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)
#* define the gr4j model
gr4j_model = HydroModel(:gr4j, components=[prod_ele, uh_flux_1, uh_flux_2, rst_ele])
# load data
df = DataFrame(CSV.File("data/gr4j/sample.csv"));
for col in names(df)[3:end]
    df[ismissing.(df[:, col]), col] .= 0.0
end
prcp_vec = df[!, "prec"]
et_vec = df[!, "pet"]
qobs_vec = df[!, "qobs"]
ts = collect(1:length(qobs_vec))

# prepare args
input = (prcp=prcp_vec, ep=et_vec)
pas = ComponentVector(
    params=(x1=320.11, x2=2.42, x3=69.63, x4=1.39),
    initstates=(soilwater=235.97, routingstore=45.47)
)
timeidx = collect(1:length(prcp_vec))
solver = HydroModels.ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3, saveat=timeidx)
result = gr4j_model(input, pas, timeidx=timeidx, solver=solver);
HydroModels.get_output_names(gr4j_model)
#! set the tunable parameters and constant parameters
tunable_pas = ComponentVector(params=(x1=320.11, x2=2.42, x3=69.63, x4=1.39))
const_pas = ComponentVector(initstates=(soilwater=235.97, routingstore=45.47))
total_pas = ComponentVector(params=(x1=320.11, x2=2.42, x3=69.63, x4=1.39), initstates=(soilwater=235.97, routingstore=45.47))
#! set the tunable parameters boundary
tunable_param_lb = [0.0, -10.0, 1.0, 0.5]
tunable_param_ub = [2000.0, 10.0, 100.0, 15.0]
#! prepare flow
output = (flow=qobs_vec,)
#! model calibration
best_pas = HydroModels.param_box_optim(
    gr4j_model,
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    input=[input],
    target=[output],
    timeidx=[timeidx],
    lb=tunable_param_lb,
    ub=tunable_param_ub,
    solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=100,
    loss_func=HydroErr.mse,
    solver=solver
)
#! use the optimized parameters for model simulation
result_opt = gr4j_model(input, HydroModels.update_ca(total_pas, ComponentVector(best_pas, getaxes(tunable_pas))), timeidx=timeidx, solver=solver)


HydroErr.mse(result_opt.flow, df[!, "qobs"])
HydroErr.mse(result.flow, df[!, "qobs"])

