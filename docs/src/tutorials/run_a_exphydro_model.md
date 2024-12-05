# Run A ExpHydro Model

## Overview

本文档将介绍如何使用已构建的exphdyro模型用于集总式(Single HRU)的计算和分布式(Multiple HURs)的计算

## Dependencies

First, lets import the required packages.

```julia
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using Plots
```

## Model Setup

接着需要对模型的参数和初始状态进行设置

### Parameter Configuration

Exphydro模型包含以下6个参数,我们需要使用ComponentVector来定义这些参数

| Parameter | Value        | Description                   |
| --------- | ------------ | ----------------------------- |
| `f`       | 0.01674478   | Infiltration parameter        |
| `Smax`    | 1709.461015  | Maximum soil water storage    |
| `Qmax`    | 18.46996175  | Maximum discharge             |
| `Df`      | 2.674548848  | Degree-day factor             |
| `Tmax`    | 0.175739196  | Maximum temperature threshold |
| `Tmin`    | -2.092959084 | Minimum temperature threshold |

```julia
params = ComponentVector(f=0.01674478, Smax=1709.461015, Qmax=18.46996175,
 Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
```

### Initial States

然后就是定义模型的初始状态,针对exphydro模型的两个计算模块分别设置其对应的状态变量`snowpack`和`soilwater`初始值

```julia
inistates = ComponentVector(snowpack=0.0, soilwater=1303.004248)
```

## Data Preparation

The test uses hydrometeorological data from "data/exphydro/01013500.csv":

```julia
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

# Run model with performance benchmarking
result = model(input, pas, config=(solver=solver, timeidx=ts), convert_to_ntp=true)

# Visualize results
plot(result.flow)
plot!(df[ts, "flow(mm)"])
```

### Performance Benchmarking For Single Node (by using BenchmarkTools.jl)

```julia
@btime model(input, pas, config=(solver=solver, timeidx=ts), convert_to_ntp=true);
```

### Multi-Node Testing (Optional)

The code includes commented sections for multi-node testing:

```julia
# Setup multiple nodes
node_num = 10
inputs = repeat([input], node_num)
ptypes = [Symbol(:node, i) for i in 1:node_num]

# Configure parameters and states for multiple nodes
params_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([params], node_num)))
init_states_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([init_states], node_num)))
pas_multi = ComponentVector(params=params_multi, initstates=init_states_multi)

# Run multi-node simulation
results = model(inputs, pas_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true)
```

### Performance Benchmarking For Multi-Node (by using BenchmarkTools.jl)

```julia
@btime model(inputs, pas_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true);
```