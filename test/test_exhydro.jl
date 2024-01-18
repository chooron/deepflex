# import lib
using CSV
using DataFrames;

# test exphydro model
include("../src/DeepFlex.jl")

# build model
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

parameters = Dict(:f => f, :Smax => Smax, :Qmax => Qmax, :Df => Df, :Tmax => Tmax, :Tmin => Tmin)
init_states = Dict(:SnowWater => 0.0, :SoilWater => 1303.004248)
ir = DeepFlex.InterceptionFilter(id="ir", parameters=parameters)
model = DeepFlex.ExpHydro(id="exp-hydro", parameters=parameters, init_states=init_states)

# load data
file_path = "data/camels/01013500.csv"

data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:100, "dayl(day)"]
prcp_vec = df[1:100, "prcp(mm/day)"]
temp_vec = df[1:100, "tmean(C)"]
flow_vec = df[1:100, "flow(mm)"]

inputs = Dict(:Prcp => prcp_vec, :Lday => lday_vec, :Temp => temp_vec)
result = DeepFlex.get_output(model, inputs)