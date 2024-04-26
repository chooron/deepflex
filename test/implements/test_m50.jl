# import lib
using CSV
using DataFrames
using BenchmarkTools

# test exphydro model
include("../../src/DeepFlex.jl")

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
time = 1:10000
lday_vec = df[time, "dayl(day)"]
prcp_vec = df[time, "prcp(mm/day)"]
temp_vec = df[time, "tmean(C)"]
flow_vec = df[time, "flow(mm)"]

# build model
model = DeepFlex.M50(name=:m50)


# pretrain model
# DeepFlex.nn_param_optim(model)


inputs = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec, time=1:1:length(lday_vec))

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = (snowwater=0.0, soilwater=1303.004248)


result = model(input, params, init_states);

result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, time, flow_vec, color=:red)
lines!(ax, time, result_df[!, :flow], color=:blue)
fig
