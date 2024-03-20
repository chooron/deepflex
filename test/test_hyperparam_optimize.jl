# import lib
using CSV
using DataFrames
using CairoMakie

# test exphydro model
include("../src/DeepFlex.jl")

# predefine the parameters
init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
paraminfos = [
    DeepFlex.BoundaryParamInfo(:SnowWater, 0.0, lb=0.0, ub=100.0),
    DeepFlex.BoundaryParamInfo(:SoilWater, 100.0, lb=0.0, ub=2000.0),
    DeepFlex.BoundaryParamInfo(:f, 0.01, lb=0.0, ub=0.1),
    DeepFlex.BoundaryParamInfo(:Smax, 100.0, lb=100.0, ub=1500.0),
    DeepFlex.BoundaryParamInfo(:Qmax, 20.0, lb=10.0, ub=50.0),
    DeepFlex.BoundaryParamInfo(:Df, 1.0, lb=0.0, ub=5.0),
    DeepFlex.BoundaryParamInfo(:Tmax, 1.0, lb=0.0, ub=3.0),
    DeepFlex.BoundaryParamInfo(:Tmin, -1.0, lb=-3.0, ub=0.0),
]
model = DeepFlex.ExpHydro(name=:exphydro)

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:10000, "dayl(day)"]
prcp_vec = df[1:10000, "prcp(mm/day)"]
temp_vec = df[1:10000, "tmean(C)"]
flow_vec = df[1:10000, "flow(mm)"]

# parameters optimization
input = ComponentVector(Prcp=prcp_vec, Lday=lday_vec, Temp=temp_vec, time=1:1:length(lday_vec))
output = ComponentVector(Flow=flow_vec)
best_params = DeepFlex.hyperparams_optimize(
    model,
    paraminfos=paraminfos,
    input=input,
    target=output)

# model = DeepFlex.ExpHydro(id="exp-hydro",
#     parameters=Dict(p.name => p.value for p in best_params),
#     initstates=Dict(p.name => p.value for p in best_params))
# result = DeepFlex.get_output(model, input=input)
# result_df = DataFrame(result)

# # plot result
# fig = Figure(size=(400, 300))
# ax = Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
# x = range(1, length(flow_vec), length=length(flow_vec))
# lines!(ax, x, flow_vec, color=:red)
# lines!(ax, x, result[:Q], color=:blue)
# fig