# HydroModels Bucket Model Implementation Example

## Overview
This document demonstrates the implementation and usage of a hydrological bucket model using the HydroModels framework. The code showcases both single-node and multi-node simulations using the ExpHydro model structure.

## Dependencies
```julia
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using HydroModels
using ModelingToolkit
```

## Model Configuration

### Parameter Setup
The model uses the following parameters for the hydrological simulation:
- `f`: 0.01674478 (Infiltration parameter)
- `Smax`: 1709.461015 (Maximum soil water storage)
- `Qmax`: 18.46996175 (Maximum discharge)
- `Df`: 2.674548848 (Degree-day factor)
- `Tmax`: 0.175739196 (Maximum temperature threshold)
- `Tmin`: -2.092959084 (Minimum temperature threshold)

### Initial States
Initial conditions for the model:
- `snowpack`: 0.0
- `soilwater`: 1303.004248

## Data Input
The model uses time series data from a CSV file located at "data/exphydro/01013500.csv" containing:
- Day length (dayl)
- Mean temperature (tmean)
- Precipitation (prcp)

## Implementation Examples

### 1. Single Node Simulation
```julia
# Setup input data
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
solver = HydroModels.ManualSolver{true}()
config = (solver=solver,)

# Convert input to required format
input_arr = Matrix(reduce(hcat, collect(input[ele.meta.inputs]))')

# Run simulation
results = ele(input_arr, pas, config=config, convert_to_ntp=true)
```

### 2. Multi-Node Simulation
```julia
# Setup for multiple nodes
node_num = 10
node_names = [Symbol(:node, i) for i in 1:node_num]

# Create parameter and state vectors for all nodes
node_params = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([params], length(node_names))))
node_initstates = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([init_states], length(node_names))))
node_pas = ComponentVector(params=node_params, initstates=node_initstates)

# Prepare input data for multiple nodes
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
node_input = permutedims(node_input, (2, 3, 1))

# Run simulation with multiple nodes
run_kwgs = (ptypes=node_names, timeidx=ts)
result = ele(node_input, node_pas, kwargs=run_kwgs)
```

## Model Structure
The model uses a bucket structure defined in exphydro.jl with two main components:
1. Surface water component (bucket_1):
   - Handles precipitation partitioning (rainfall/snowfall)
   - Computes potential evapotranspiration
   - Manages snowmelt processes

2. Soil water component:
   - Manages soil water storage
   - Computes actual evaporation
   - Generates baseflow and surface flow

## Usage Notes
1. The code demonstrates flexibility in handling both single-node and multi-node simulations
2. Input data should be properly formatted with required columns (dayl, tmean, prcp)
3. Parameters and initial states can be adjusted based on specific catchment characteristics
4. The model uses a manual solver for time-stepping

## Time Series Processing
The example processes 10,000 time steps (ts = collect(1:10000)) and can be adjusted based on data availability and simulation requirements.