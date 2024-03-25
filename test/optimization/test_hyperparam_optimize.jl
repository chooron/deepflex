# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using OptimizationOptimisers

# test exphydro model
include("../../src/DeepFlex.jl")

# predefine the parameters
init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
search_params = [
    # DeepFlex.BoundaryParamInfo(:snowwater, 0.01, lb=0.0, ub=10.0),
    # DeepFlex.BoundaryParamInfo(:soilwater, 1000.01, lb=100.0, ub=1500.0),
    DeepFlex.BoundaryParamInfo(:f, 0.01, lb=0.0, ub=0.1),
    DeepFlex.BoundaryParamInfo(:Smax, 100.0, lb=100.0, ub=1500.0),
    DeepFlex.BoundaryParamInfo(:Qmax, 20.0, lb=10.0, ub=50.0),
    DeepFlex.BoundaryParamInfo(:Df, 1.0, lb=0.0, ub=5.0),
    DeepFlex.BoundaryParamInfo(:Tmax, 1.0, lb=0.0, ub=3.0),
    DeepFlex.BoundaryParamInfo(:Tmin, -1.0, lb=-3.0, ub=0.0),
]
const_params = (snowwater=0.0,soilwater=1300.0)
# const_params = NamedTuple()
model = DeepFlex.ExpHydro(name=:exphydro)

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec, time=1:1:length(lday_vec))
output = (flow=flow_vec,)

best_params = DeepFlex.hyperparams_optimize(
    model,
    search_params=search_params,
    const_params=const_params,
    input=input,
    target=output,
    # solve_alg=OptimizationOptimisers.Adam(0.01)
)

total_params = merge(best_params, const_params)
result = DeepFlex.get_output(model, input=input,
    parameters=total_params[model.param_names],
    init_states=total_params[model.state_names])
