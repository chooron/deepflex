# Running a Single Component in HydroModels

This tutorial will guide you through the process of running a single component in the HydroModels package. We'll use the `SurfaceStorage` component from the `ExpHydro` model as an example.

## Step 1: Import Required Modules

First, let's import the necessary modules:

```julia
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using ModelingToolkit
```

## Step 2: Load Data

```julia
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
```

## Step 3: Define Parameters and Initial States

```julia
# load predefined parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=init_states)
```

## Step 4: Define the Component

```julia
# load predefined component in implemented models
ele = HydroModels.ExpHydro.SurfaceStorage(name=:sf)
```

## Step 5: Run the Component

```julia
bucket_input_names = HydroModels.get_input_names(ele)
# we need to sort the input variables to make sure they are in the correct order.
# the input of the component is a matrix with the shape of (var nm * ts len).
input_matrix = Matrix(reduce(hcat, collect(input[bucket_input_names]))') # (var nm * ts len)
# define solver
solver = HydroModels.ODESolver()
results = ele(input_matrix, pas, timeidx=ts, solver=solver)
```

## Step 6: Show Results

the output results include both the outputs and states of the component.

```julia
bucket_output_names = HydroModels.get_output_names(ele)
bucket_state_names = HydroModels.get_state_names(ele)
results_ntp = NamedTuple{Tuple(vcat(bucket_output_names, bucket_state_names))}(eachslice(results, dims=1))
results_df = DataFrame(results_ntp)
# show data
@info first(results_df, 5)
```
