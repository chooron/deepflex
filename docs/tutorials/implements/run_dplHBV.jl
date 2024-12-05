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
using DataInterpolations
using Optimization
using OptimizationOptimisers
using SciMLSensitivity
using JLD2
using Plots
include("../../../src/HydroModels.jl")
include("../models/dplHBV.jl")
include("loss_functions.jl")

#* load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
ts = collect(1:10000)
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
dayl_vec = df[ts, "dayl(day)"]
pet_vec = @. 29.8 * dayl_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)
qobs_vec = df[ts, "flow(mm)"]
input = (prcp=prcp_vec, pet=pet_vec, temp=temp_vec)
#* prepare parameters
psnn_ps, psnn_st = Lux.setup(StableRNG(123), params_nn)
psnn_ps_ca = ComponentVector(psnn_ps)
psnn_ps_vec = Vector(ComponentVector(psnn_ps))
params = ComponentVector(TT=0.0, CFMAX=5.0, CWH=0.1, CFR=0.05, FC=200.0, LP=0.6, k0=0.06, k1=0.2, k2=0.1, PPERC=2, UZL=10)
nns = ComponentVector(NamedTuple{Tuple([nn_wrapper.meta.name])}([psnn_ps_vec]))
init_states = ComponentVector(suz=0.0, slz=0.0, soilwater=0.0, meltwater=0.0, snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states, nn=nns)
#* define config
config = (solver=HydroModels.ODESolver(), timeidx=ts, interp=LinearInterpolation)
#* run model
# @btime result = model(input, pas, timeidx=ts, convert_to_ntp=true)

# sum((result.q .- qobs_vec) .^ 2)  / sum((qobs_vec .- mean(qobs_vec)) .^2) 
# # # #! set the tunable parameters boundary
# # #! prepare flow
# tunable_pas = ComponentVector(params=params, nn=(pnn=psnn_ps_vec,))
# const_pas = ComponentVector(initstates=init_states)

# model_grad_opt = HydroModels.GradOptimizer(component=model, maxiters=100, adtype=AutoZygote(), solve_alg=Adam(1e-2))
# config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autodiff=true)), timeidx=ts, interp=LinearInterpolation)
# hbv_hydro_opt_params, loss_df = model_grad_opt(
#     [input], [(q=qobs_vec,)],
#     tunable_pas=tunable_pas,
#     const_pas=const_pas,
#     return_loss_df=true
# )

# mse_value, fhv_value, nse_value = mse(qobs_vec, result.q), fhv(qobs_vec, result.q, 0.1), nse(qobs_vec, result.q)

# plot(result.q, label="simulated-before-optimal")
# plot!(re_result.q, label="simulated-after-optimal")
# plot!(qobs_vec, label="observed")

# hbv_hydro_opt_params = load("doc/tutorials/implements/save/dplHBV_opt.jld2", "opt_params")
# re_result = model(input, hbv_hydro_opt_params, timeidx=ts, convert_to_ntp=true)
# plot(re_result.q, label="simulated")
# plot!(qobs_vec, label="observed")
# 1 - sum((re_result.q .- qobs_vec) .^ 2)  / sum((qobs_vec .- mean(qobs_vec)) .^2)

# # save result
# save("doc/tutorials/implements/save/dplHBV_opt.jld2", "loss_df", loss_df, "opt_params", hbv_hydro_opt_params)