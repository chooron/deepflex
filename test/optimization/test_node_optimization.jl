# import lib
using CSV
using DataFrames
# using CairoMakie
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools
using NamedTupleTools

# test exphydro model
include("../../src/DeepFlex.jl")

# predefine the parameters
# init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
init_param_axes = Axis(:f, :Smax, :Qmax, :Df, :Tmax, :Tmin)
init_params_list = [0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084]

default_params = ComponentVector(exphydro=(unit=ComponentVector(init_params_list, init_param_axes), route=ComponentVector(), weight=1.0))

p1 = ComponentVector(ep1=(unit=1.0,),) # ComponentVector
p2 = ComponentVector(ep2=(unit=1.0,),) # ComponentVector
p3 = ComponentVector(ep3=(weight=1.0,),) # ComponentVector

ComponentVector(p1[:ep]; p2[:ep]...)
ComponentVector(p1; p3...)
lb_list = [0.0, 100.0, 10.0, 0.0, 0.0, -3.0]
ub_list = [0.1, 1500.0, 50.0, 5.0, 3.0, 0.0]

tunable_param_values = ComponentVector(init_params, tunable_param_axes)
tunable_param_lb = ComponentVector(lb_list, tunable_param_axes)
tunable_param_ub = ComponentVector(ub_list, tunable_param_axes)
const_params_values = ComponentVector(snowwater=0.0, soilwater=1300.0)

model = DeepFlex.ExpHydro_Node(name=:exphydro)

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

best_params = DeepFlex.param_box_optim(
    model,
    tunable_params=tunable_param_values,
    const_params=const_params_values,
    input=input,
    target=output,
    lb=tunable_param_lb,
    ub=tunable_param_ub,
)

total_params = merge(best_params, const_params)
reulst = model(input, total_params[model.nameinfo.param_names], total_params[model.nameinfo.state_names])

