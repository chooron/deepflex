# import lib
using CSV
using DataFrames
using CairoMakie
using BenchmarkTools
using ComponentArrays
# using LumpedHydro
# test exphydro model
include("../../src/LumpedHydro.jl")
# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
flow_vec = df[ts, "flow(mm)"]

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
unit_params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_init_states = (snowwater=0.0, soilwater=1303.004248)

pas = ComponentVector((params=unit_params, initstates=unit_init_states, weight=1.0))

model = LumpedHydro.ExpHydro.Unit(name=:exphydro, mtk=false)

input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
solver = LumpedHydro.ODESolver(reltol=1e-3, abstol=1e-3)
# solver = LumpedHydro.DiscreteSolver()
result = model(input, pas, timeidx=ts, solver=solver);
# result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts, flow_vec, color=:red)
lines!(ax, ts, result.flow, color=:blue)

fig