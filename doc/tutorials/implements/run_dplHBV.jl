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
using OptimizationOptimisers
using HydroErrors
include("../../../src/HydroModels.jl")
include("../models/dplHBV.jl")

# load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
ts = collect(1:100)
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
dayl_vec = df[ts, "dayl(day)"]
qobs_vec = df[ts, "flow(mm)"]
psnn_ps, psnn_st = Lux.setup(StableRNG(123), params_nn)
psnn_ps_ca = ComponentVector(psnn_ps)
psnn_ps_vec = Vector(ComponentVector(psnn_ps))
pet_vec = @. 29.8 * dayl_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)
params = ComponentVector(TT=0.0, CFMAX=5.0, CWH=0.1, CFR=0.05, FC=200.0, LP=0.6, k0=0.06, k1=0.2, k2=0.1, PPERC=2, UZL=10)
# nns = ComponentVector(NamedTuple{Tuple([params_nn.name])}([psnn_ps_vec]))
init_states = ComponentVector(suz=0.0, slz=0.0, soilwater=0.0, meltwater=0.0, snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states, nn=(pnn=psnn_ps_vec,))
input = (prcp=prcp_vec, pet=pet_vec, temp=temp_vec)
result = model(input, pas, timeidx=ts, convert_to_ntp=true)


# #! set the tunable parameters boundary
# #! prepare flow
tunable_pas = ComponentVector(params=params, nn=(pnn=psnn_ps_vec,))
const_pas = ComponentVector(initstates=init_states)

model_grad_opt = HydroModels.GradOptimizer(component=model, maxiters=100,adtype=AutoZygote(), solve_alg=Adam(1e-2))
hbv_hydro_opt_params, loss_df = model_grad_opt(
    [input], [(q=qobs_vec,)],
    tunable_pas=tunable_pas,
    const_pas=const_pas,
    return_loss_df=true
)
# plot(loss_df[!,:loss])