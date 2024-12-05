using CSV
using Lux
using LuxCore
using Random
using DataFrames
using Symbolics
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using BenchmarkTools
using StableRNGs
using Optimization
using OptimizationBBO
using HydroErrors
include("../../../src/HydroModels.jl")
include("../models/HBV.jl")

# load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
ts = collect(1:10000)
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
dayl_vec = df[ts, "dayl(day)"]
qobs_vec = df[ts, "flow(mm)"]
pet_vec = @. 29.8 * dayl_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)
params = ComponentVector(TT=0.0, CFMAX=5.0, CWH=0.1, CFR=0.05, FC=200.0, LP=0.6, BETA=3.0, k0=0.06, k1=0.2, k2=0.1, PPERC=2, UZL=10)
init_states = ComponentVector(suz=0.0, slz=0.0, soilwater=0.0, meltwater=0.0, snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states)
input = (prcp=prcp_vec, pet=pet_vec, temp=temp_vec)
result = model(input, pas, timeidx=ts, convert_to_ntp=true)
#! set the tunable parameters boundary
lower_bounds = [-1.5, 1, 0.0, 0.0, 50.0, 0.3, 1.0, 0.05, 0.01, 0.001, 0.0, 0.0]
upper_bounds = [1.2, 8.0, 0.2, 0.1, 500.0, 1.0, 6.0, 0.5, 0.3, 0.15, 3.0, 70.0]
#! prepare flow
tunable_pas = ComponentVector(params=params)
const_pas = ComponentVector(initstates=init_states)

model_hydro_opt = HydroModels.HydroOptimizer(component=model, maxiters=100)
hbv_hydro_opt_params, loss_df = model_hydro_opt(
    [input], [(q=qobs_vec,)],
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    lb=lower_bounds,
    ub=upper_bounds,
    return_loss_df=true
)
plot(loss_df[!,:loss])

# #* ComponentVector{Float64}(params = (TT = -1.2223657527438707, CFMAX = 2.201359793941345, CWH = 0.022749518921432663, CFR = 0.058335602629828544, FC = 160.01327559173077, LP = 0.7042581781418978, 
# #* beta = 5.580695551758287, k0 = 0.0500023960318018, k1 = 0.04573064980956475, k2 = 0.14881856483902567, PERC = 1.3367222956722589, UZL = 44.059927907190016))
# #! model calibration
# best_pas = HydroModels.param_box_optim(
#     hbv_model,
#     tunable_pas=ComponentVector(params=params),
#     const_pas=ComponentVector(initstates=init_states),
#     input=input,
#     target=output,
#     timeidx=ts,
#     lb=lower_bounds,
#     ub=upper_bounds,
#     solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
#     maxiters=1000,
# )

# result = hbv_model(input, HydroModels.merge_ca(pas, best_pas), timeidx=ts)
# HydroModels.nse(result.flow, qobs_vec)