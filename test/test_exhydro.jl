# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays

# test exphydro model
include("../src/DeepFlex.jl")

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

parameters = ComponentVector(; Dict(:f => f, :Smax => Smax, :Qmax => Qmax, :Df => Df, :Tmax => Tmax, :Tmin => Tmin)...)
init_states = ComponentVector(; Dict(:SnowWater => 0.0, :SoilWater => 1303.004248)...)
model = DeepFlex.ExpHydro(id="exp-hydro", parameters=parameters, init_states=init_states)

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec)
result = DeepFlex.get_output(model, input=inputs, step=false)

# result_df = DataFrame(result)

# # plot result
# fig = Figure(size=(400, 300))
# ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
# x = range(1, 1000, length=1000)
# lines!(ax, x, flow_vec, color=:red)
# lines!(ax, x, result[:Q], color=:blue)
# fig