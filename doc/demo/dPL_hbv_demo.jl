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
include("../../src/HydroModels.jl")

# load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
ts = collect(1:10000)
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
dayl_vec = df[ts, "dayl(day)"]
qobs_vec = df[ts, "flow(mm)"]

#! parameters in the HBV-light model
@parameters TT CFMAX CFR CWH FC beta LP PERC k0 k1 k2 UZL

#! hydrological flux in the dpl-hbv model
@variables prcp temp lday pet rainfall snowfall refreeze melt snowinfil
@variables snowpack meltwater soilwater
@variables soilwetfrac recharge excess evapfrac evap
@variables upperzone lowerzone perc q0 q1 q2 flow

#* dynamic parameters
@variables beta gamma

SimpleFlux = HydroModels.SimpleFlux
NeuralFlux = HydroModels.NeuralFlux
StateFlux = HydroModels.StateFlux
HydroModel = HydroModels.HydroModel
HydroBucket = HydroModels.HydroBucket
step_func = x -> ifelse(x > 0.0, 1.0, 0.0)

#* θ是regional data所以可以针对实测和预测结果模拟优化出来, 
#* 而γ和β是根据观测数据预测得到的随时间变化的参数,
#* 因此对于某个单个流域的解决方式就是,仅通过nn预测这两个参数,然后与θ一同进行优化
#* 如何嵌入LSTM到模型中? 模型输入是二维的矩阵而非一维向量了,但是timeidx会发生改变,
#* 先考虑用mlp这种简单的比较好

#! define the snow pack reservoir
snow_funcs = [
    SimpleFlux([prcp, temp] => [snowfall, rainfall], [TT], exprs=[step_func(TT - temp) * prcp, step_func(temp - TT) * prcp]),
    SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    SimpleFlux([snowpack, temp] => [melt], [TT, CFMAX], exprs=[min(snowpack, CFMAX * max(0.0, temp - TT))]),
    SimpleFlux([meltwater, temp] => [refreeze], [TT, CFMAX, CFR], exprs=[min(meltwater, CFR * CFMAX * max(0.0, TT - temp))]),
    SimpleFlux([meltwater, snowpack] => [snowinfil], [CWH], exprs=[max(0.0, meltwater - CWH * snowpack)]),
]
snow_dfuncs = [StateFlux([snowfall, refreeze] => [melt], snowpack), StateFlux([melt] => [refreeze, snowinfil], meltwater)]
snow_ele = HydroBucket(name=:hbv_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

thetas_nn = Lux.Chain(
    Lux.Dense(4 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 2, Lux.leakyrelu),
    name=:thetas
)
recharge_nn = Lux.Chain(
    Lux.Dense(2 => 16, Lux.tanh),
    Lux.Dense(16 => 16, Lux.leakyrelu),
    Lux.Dense(16 => 1, Lux.leakyrelu),
    name=:recharge
)

thetas_nn_ps = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(123), thetas_nn)))
recharge_nn_ps = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(123), recharge_nn)))

#! define the soil water reservoir
soil_funcs = [
    NeuralFlux([prcp, temp, pet, soilwater] => [beta, gamma], thetas_nn),
    SimpleFlux([soilwater, beta] => [soilwetfrac], [FC], exprs=[clamp(abs(soilwater / FC)^beta, 0.0, 1.0)]),
    NeuralFlux([rainfall, snowinfil, soilwetfrac] => [recharge], recharge_nn),
    SimpleFlux([soilwater] => [excess], [FC], exprs=[max(0.0, soilwater - FC)]),
    SimpleFlux([soilwater] => [evapfrac], [LP, FC], exprs=[soilwater / (LP * FC)]),
    SimpleFlux([soilwater, pet, evapfrac, gamma] => [evap], exprs=[min(soilwater, pet * (abs(evapfrac)^gamma))]),
]

soil_dfuncs = [StateFlux([rainfall, snowinfil] => [evap, excess, recharge], soilwater)]
soil_ele = HydroBucket(name=:hbv_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

#! define the upper and lower subsurface zone 
zone_funcs = [
    SimpleFlux([soilwater, upperzone, evapfrac] => [perc], [PERC], exprs=[min(soilwater, PERC)]),
    SimpleFlux([upperzone] => [q0], [k0, UZL], exprs=[max(0.0, (upperzone - UZL) * k0)]),
    SimpleFlux([upperzone] => [q1], [k1], exprs=[upperzone * k1]),
    SimpleFlux([lowerzone] => [q2], [k2], exprs=[lowerzone * k2]),
    SimpleFlux([q0, q1, q2] => [flow], exprs=[q0 + q1 + q2]),
]

zone_dfuncs = [StateFlux([recharge, excess] => [perc, q0, q1], upperzone), StateFlux([perc] => [q2], lowerzone)]
zone_ele = HydroBucket(name=:hbv_zone, funcs=zone_funcs, dfuncs=zone_dfuncs)

#! define the Exp-Hydro model
hbv_model = HydroModel(name=:hbv, components=[snow_ele, soil_ele, zone_ele]);

params = ComponentVector(TT=0.0, CFMAX=5.0, CWH=0.1, CFR=0.05, FC=200.0, LP=0.6, beta=3.0, k0=0.06, k1=0.2, k2=0.1, PERC=2, UZL=10)
init_states = ComponentVector(upperzone=0.0, lowerzone=0.0, soilwater=0.0, meltwater=0.0, snowpack=0.0)
nn_params = ComponentVector(thetas=thetas_nn_ps, recharge=recharge_nn_ps)
pas = ComponentVector(params=params, initstates=init_states, nn=nn_params)
input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
config = (solver=HydroModels.ODESolver(), timeidx=ts)
# @btime result = hbv_model(input, pas, config=config)

#! prepare flow
output = (flow=qobs_vec,)
#! model calibration
#* ComponentVector{Float64}(params = (TT = -1.2223657527438707, CFMAX = 2.201359793941345, CWH = 0.022749518921432663, CFR = 0.058335602629828544, FC = 160.01327559173077, LP = 0.7042581781418978, 
#* beta = 5.580695551758287, k0 = 0.0500023960318018, k1 = 0.04573064980956475, k2 = 0.14881856483902567, PERC = 1.3367222956722589, UZL = 44.059927907190016))
best_pas = HydroModels.param_grad_optim(
    hbv_model,
    tunable_pas=ComponentVector(params=params, nn=nn_params),
    const_pas=ComponentVector(initstates=init_states),
    input=[input],
    target=[output],
    maxiters=100,
    loss_func=HydroErrors.mse,
)