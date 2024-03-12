# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using BenchmarkTools

# test exphydro model
include("../src/DeepFlex.jl")

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

# build solver
solver = DeepFlex.ODESolver(dt=1, saveat=1:1:length(lday_vec))

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
parameters = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(SnowWater=0.0, SoilWater=1303.004248)

model = DeepFlex.ExpHydro(name="exp-hydro")

DeepFlex.set_parameters!(model, parameters)
DeepFlex.set_init_states!(model, init_states)
DeepFlex.set_solver!(model, solver)

inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec)
result = DeepFlex.get_output(model, input=inputs, step=true, sensealg=DeepFlex.default_ode_sensealg);

result_df = DataFrame(Dict(k => result[k] for k in keys(result)))
result_df[!, :snowwater] = model.elements[1].states[:SnowWater]
result_df[!, :soilwater] = model.elements[2].states[:SoilWater]

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 1000, length=1000)
lines!(ax, x, flow_vec, color=:red)
lines!(ax, x, result_df[!, :Flow], color=:blue)
fig