using CSV            
using DataFrames     
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using Plots
using HydroModels
using HydroModelTools
using OrdinaryDiffEq
include("../models/ExpHydro.jl")
println(exphydro_model)

# Define model parameters
params = ComponentVector(
    f = 0.01674478,    # Infiltration parameter
    Smax = 1709.461015, # Maximum soil water storage
    Qmax = 18.46996175, # Maximum discharge
    Df = 2.674548848,   # Degree-day factor
    Tmax = 0.175739196, # Maximum temperature threshold
    Tmin = -2.092959084 # Minimum temperature threshold
)

# Define initial states
inistates = ComponentVector(
    snowpack = 0.0,           # Snow accumulation
    soilwater = 1303.004248   # Soil water content
)

pas = ComponentVector(params=params,initstates=inistates)

# Read input data
file_path = "../data/exphydro/01013500.csv"
df = DataFrame(CSV.File(file_path))
ts = collect(1:10000)  # Time series length
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(exphydro_model)]))')

# Configure solver
# solver = HydroModels.ManualSolver{true}()
solver = ODESolver()

# Run model with benchmarking
@btime result = exphydro_model(input, pas, config=(solver=solver, timeidx=ts));

@btime result1 = exphydro_model(input, pas, config=(solver=ODESolver(alg=Tsit5(), abstol=1e-3, reltol=1e-3), timeidx=ts));
@btime result2 = exphydro_model(input, pas, config=(solver=DiscreteSolver(), timeidx=ts));

# convert to NamedTuple
model_total_output_names = vcat(HydroModels.get_state_names(exphydro_model), HydroModels.get_output_names(exphydro_model))
result_ntp = NamedTuple{Tuple(model_total_output_names)}(eachslice(result, dims=1))

DataFrame(result_ntp)
# Visualize results
plot(result_ntp.flow, label="Simulated")
plot!(df[ts, "flow(mm)"], label="Observed")

savefig("exphydro_predict.png")

# Setup multiple nodes
node_num = 10  # Number of nodes
inputs = repeat(reshape(input, size(input)[1], 1, size(input)[2]), 1, node_num, 1)
ptypes = [Symbol(:node, i) for i in 1:node_num]

# Configure parameters and states for multiple nodes
params_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([params], node_num)))
init_states_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([inistates], node_num)))
pas_multi = ComponentVector(params=params_multi, initstates=init_states_multi)
# Run multi-node simulation
results = exphydro_model(inputs, pas_multi, config=(solver=solver, timeidx=ts))