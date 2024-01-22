# import lib
using CSV
using DataFrames

# test exphydro model
include("../src/DeepFlex.jl")

# predefine the parameters
paraminfos = [
    DeepFlex.ParamInfo(:SnowWater, 0.0, lb=0.0, ub=100.0),
    DeepFlex.ParamInfo(:SoilWater, 100.0, lb=0.0, ub=2000.0),
    DeepFlex.ParamInfo(:f, 0.01, lb=0.0, ub=0.1),
    DeepFlex.ParamInfo(:Smax, 100, lb=100, ub=1500),
    DeepFlex.ParamInfo(:Qmax, 20, lb=10, ub=50),
    DeepFlex.ParamInfo(:Df, 1.0, lb=0.0, ub=5.0),
    DeepFlex.ParamInfo(:Tmax, 1.0, lb=0.0, ub=3.0),
    DeepFlex.ParamInfo(:Tmin, -1.0, lb=-3.0, ub=0.0),
]
model = DeepFlex.ExpHydro(id="exp-hydro", paraminfos=paraminfos)

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

# parameters optimization
input = Dict(:Prcp => prcp_vec, :Lday => lday_vec, :Temp => temp_vec)
output = Dict(:Q => flow_vec)
best_params = DeepFlex.hyper_params_optimize(model, paraminfos, input, output)


# # plot result
# fig = Figure(size=(400, 300))
# ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
# x = range(1, 1000, length=1000)
# lines!(ax, x, flow_vec, color=:red)
# lines!(ax, x, result[:Q], color=:blue)
