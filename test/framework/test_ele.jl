# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools

include("../../src/DeepFlex.jl")

ele = DeepFlex.ExpHydro.Surface(name=:sf, mtk=true)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:100
input = ComponentVector(time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
solver = DeepFlex.ODESolver()
# sys = DeepFlex.setup_input(ele, input=input, time=ts)
# prob = ODEProblem(sys, [], (1, 10000), [])
# new_p = SciMLStructures.replace(Tunable(), prob.p, [-2.09, 2.6745, 0.1757])

# new_prob = remake(prob,
#     p=[sys.sf_surf_base_sys.Df => Df, sys.sf_surf_base_sys.Tmax => Tmax, sys.sf_surf_base_sys.Tmin => Tmin],
#     u0=[sys.sf_surf_base_sys.snowwater => 0.0])
# ModelingToolkit.MTKParameters
results = ele(input, pas, solver=solver)
vcat(results.u...)

# DeepFlex.@simpleflux([:temp, :lday], "pet", Symbol[])