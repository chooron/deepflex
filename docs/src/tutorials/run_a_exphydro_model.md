# ExpHydro Model Tutorial

## Overview

This tutorial demonstrates how to use the ExpHydro model for hydrological calculations, including:
- Lumped calculation (Single HRU)
- Distributed calculation (Multiple HRUs)

## Dependencies

First, import the required Julia packages:

```julia
using CSV            # Data reading
using DataFrames     # Data processing
using ComponentArrays # Parameter organization
using BenchmarkTools # Performance testing
using NamedTupleTools # Named tuple utilities
using Plots         # Result visualization
```

## Model Setup

### Parameter Configuration

The ExpHydro model contains 6 key parameters defined using ComponentVector:

| Parameter | Value        | Description                   |
| --------- | ------------ | ----------------------------- |
| `f`       | 0.01674478  | Infiltration parameter        |
| `Smax`    | 1709.461015 | Maximum soil water storage    |
| `Qmax`    | 18.46996175 | Maximum discharge             |
| `Df`      | 2.674548848 | Degree-day factor             |
| `Tmax`    | 0.175739196 | Maximum temperature threshold |
| `Tmin`    | -2.092959084| Minimum temperature threshold |

```julia
# Define model parameters
params = ComponentVector(
    f = 0.01674478,    # Infiltration parameter
    Smax = 1709.461015, # Maximum soil water storage
    Qmax = 18.46996175, # Maximum discharge
    Df = 2.674548848,   # Degree-day factor
    Tmax = 0.175739196, # Maximum temperature threshold
    Tmin = -2.092959084 # Minimum temperature threshold
)
```

### Initial States

Set initial states for the two calculation modules:

```julia
# Define initial states
inistates = ComponentVector(
    snowpack = 0.0,           # Snow accumulation
    soilwater = 1303.004248   # Soil water content
)
```

## Data Preparation

The model uses hydrometeorological data from "data/exphydro/01013500.csv":

```julia
# Read input data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)  # Time series length
```

Input variables include:

- Day length (dayl)
- Mean temperature (tmean)
- Precipitation (prcp)

## Model Testing

### Single Node Testing

```julia
# Configure solver
solver = HydroModels.ManualSolver()

# Run model with benchmarking
result = model(input, params, config=(solver=solver, timeidx=ts), convert_to_ntp=true)

# Visualize results
plot(result.flow, label="Simulated")
plot!(df[ts, "flow(mm)"], label="Observed")
```

### Single Node Benchmarking

```julia
# Performance testing using BenchmarkTools
@btime model(input, params, config=(solver=solver, timeidx=ts), convert_to_ntp=true);
```

### Multi-Node Testing

Multi-node setup for distributed computing:

```julia
# Setup multiple nodes
node_num = 10  # Number of nodes
inputs = repeat([input], node_num)
ptypes = [Symbol(:node, i) for i in 1:node_num]

# Configure parameters and states for multiple nodes
params_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([params], node_num)))
init_states_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([inistates], node_num)))

# Run multi-node simulation
results = model(inputs, params_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true)
```

### Multi-Node Benchmarking

```julia
# Multi-node performance testing
@btime model(inputs, params_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true);