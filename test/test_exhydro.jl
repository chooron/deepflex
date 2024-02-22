# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays

# test exphydro model
include("../src/DeepFlex.jl")

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
parameters = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(SnowWater=0.0, SoilWater=1303.004248)
model = DeepFlex.ExpHydro(name="exp-hydro", parameters=parameters, init_states=init_states)

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:10000, "dayl(day)"]
prcp_vec = df[1:10000, "prcp(mm/day)"]
temp_vec = df[1:10000, "tmean(C)"]
flow_vec = df[1:10000, "flow(mm)"]

inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec)
result = DeepFlex.get_output(model, input=inputs, step=false, sensealg=DeepFlex.default_ode_sensealg)

result_df = DataFrame(Dict(k => result[k] for k in keys(result)))

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 10000, length=10000)
lines!(ax, x, flow_vec, color=:red)
lines!(ax, x, result[:Surfaceflow], color=:blue)
lines!(ax, x, result[:Baseflow], color=:blue)
fig

# plot states
states = DeepFlex.get_states(model) # , state_names=Set([:SnowWater,:SoilWater])

result_df[!,:SoilWater] = states[:SoilWater]
result_df[!,:SnowWater] = states[:SnowWater]

CSV.write("data/cache/01013500.csv", result_df);