# 导入模块
using CSV
using DataFrames
using ComponentArrays
using StructArrays
using BenchmarkTools
using NamedTupleTools
using OrdinaryDiffEq

include("../../src/LumpedHydro.jl")

ele = LumpedHydro.ExpHydro.Surface(name=:sf, mtk=true, ptype=:discrete)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
ps = [f, Smax, Qmax, Df, Tmax, Tmin]
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

function test_build_system()
    build_sys = LumpedHydro.setup_input(ele.system, input, ts, LumpedHydro.get_input_names(ele), :test)
    init_prob = ODEProblem(build_sys, Pair[], (1, 100), [])
    init_prob
end
solver = LumpedHydro.DiscreteSolver()
results = ele(input, pas, timeidx=ts, solver=solver)