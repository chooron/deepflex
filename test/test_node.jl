# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis

# test exphydro model
include("../src/DeepFlex.jl")

# load data
file_path = "data/cache/01013500.csv"

exphydro_df = DataFrame(CSV.File("data/cache/01013500.csv"));
data_df = DataFrame(CSV.File("data/camels/01013500.csv"));

lday_vec = data_df[1:1000, "dayl(day)"]
prcp_vec = data_df[1:1000, "prcp(mm/day)"]
temp_vec = data_df[1:1000, "tmean(C)"]
flow_vec = data_df[1:1000, "flow(mm)"]

snowwater_vec = exphydro_df[1:1000, "SnowWater"]
soilwater_vec = exphydro_df[1:1000, "SoilWater"]
evap_vec = exphydro_df[1:1000, "Evap"]
pred_flow_vec = exphydro_df[1:1000, "Flow"]

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
parameters = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(SnowWater=0.0, SoilWater=1303.004248)
model = DeepFlex.M50(id="M50", parameters=parameters, init_states=init_states)

train_inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec,
    SnowWater=snowwater_vec, SoilWater=soilwater_vec, Evap=evap_vec, Flow=pred_flow_vec)
DeepFlex.pretrain!(model, input=train_inputs)
@info "pretrain complete!"

ode_inputs = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec)
result = DeepFlex.get_output(model, input=ode_inputs, step=false, sensealg=DeepFlex.default_node_sensealg, reltol=1e-3, abstol=1e-3)

result_df = DataFrame(Dict(k => result[k] for k in keys(result)))
# plot result
fig = Figure(size=(400, 300))
ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 1000, length=1000)
lines!(ax, x, flow_vec, color=:red)
lines!(ax, x, result[:Flow], color=:blue)
fig

