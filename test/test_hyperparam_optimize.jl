# import lib
using CSV
using DataFrames
using CairoMakie

# test exphydro model
include("../src/DeepFlex.jl")

# predefine the parameters
init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
paraminfos = [
    DeepFlex.ParamInfo(:SnowWater, 0.0, lb=0.0, ub=100.0),
    DeepFlex.ParamInfo(:SoilWater, 100.0, lb=0.0, ub=2000.0),
    DeepFlex.ParamInfo(:f, 0.01, lb=0.0, ub=0.1),
    DeepFlex.ParamInfo(:Smax, 100.0, lb=100.0, ub=1500.0),
    DeepFlex.ParamInfo(:Qmax, 20.0, lb=10.0, ub=50.0),
    DeepFlex.ParamInfo(:Df, 1.0, lb=0.0, ub=5.0),
    DeepFlex.ParamInfo(:Tmax, 1.0, lb=0.0, ub=3.0),
    DeepFlex.ParamInfo(:Tmin, -1.0, lb=-3.0, ub=0.0),
]
model = DeepFlex.ExpHydro(id="exp-hydro",
    parameters=Dict(:f => 0.01, :Smax => 100.0, :Qmax => 20.0, :Df => 1.0, :Tmax => 1.0, :Tmin => -1.0),
    initstates=Dict(:SnowWater => 0.0, :SoilWater => 100.0))

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[!, "dayl(day)"]
prcp_vec = df[!, "prcp(mm/day)"]
temp_vec = df[!, "tmean(C)"]
flow_vec = df[!, "flow(mm)"]

# parameters optimization
input = Dict(:Prcp => prcp_vec, :Lday => lday_vec, :Temp => temp_vec)
output = Dict(:Q => flow_vec)
best_params = DeepFlex.hyper_params_optimize(model, paraminfos, input, output)

model = DeepFlex.ExpHydro(id="exp-hydro",
    parameters=Dict(p.name=>p.value for p in best_params),
    initstates=Dict(p.name=>p.value for p in best_params))
result = DeepFlex.get_output(model, input=input)
result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, length(flow_vec), length=length(flow_vec))
lines!(ax, x, flow_vec, color=:red)
lines!(ax, x, result[:Q], color=:blue)
fig