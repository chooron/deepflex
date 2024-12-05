# Introduction to the implementation of the M50 model and its training process

## Step 1: Import packages
```julia
using CSV
using DataFrames
using Lux
using ModelingToolkit
using StableRNGs
using ComponentArrays
using DataInterpolations
using Statistics
using BenchmarkTools
using Plots
using OptimizationOptimisers
using SciMLSensitivity
using HydroModels
```

## Step 2: Define the model
```julia
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
#! parameters in the Exp-Hydro model
@parameters Tmin Tmax Df Smax f Qmax
#! parameters in normalize flux
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

#! hydrological flux in the Exp-Hydro model
@variables prcp temp lday pet rainfall snowfall
@variables snowpack soilwater lday pet
@variables melt log_evap_div_lday log_flow
@variables norm_snw norm_slw norm_temp norm_prcp

#! define the snow pack reservoir
snow_funcs = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the ET NN and Q NN
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)
ep_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), ep_nn))))
q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
q_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))))

#! get init parameters for each NN
ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

#! define the soil water reservoir
soil_funcs = [
    #* normalize
    HydroFlux([snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [snowpack_mean, soilwater_mean, prcp_mean, temp_mean, snowpack_std, soilwater_std, prcp_std, temp_std],
        exprs=[(var - mean) / std for (var, mean, std) in zip([snowpack, soilwater, prcp, temp],
            [snowpack_mean, soilwater_mean, prcp_mean, temp_mean],
            [snowpack_std, soilwater_std, prcp_std, temp_std]
        )]),
    ep_nn_flux,
    q_nn_flux,
]

state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr)]
soil_ele = HydroBucket(name=:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)
#! define the Exp-Hydro model
m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele]);
```

## Step 3: Load data
```julia
# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)
lday_vec = df[ts, "Lday"]
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
flow_vec = df[ts, "Flow"]

log_flow_vec = log.(flow_vec)
log_evap_div_lday_vec = log.(df[ts, "Evap"] ./ lday_vec)
norm_prcp_vec = (prcp_vec .- mean(prcp_vec)) ./ std(prcp_vec)
norm_temp_vec = (temp_vec .- mean(temp_vec)) ./ std(temp_vec)
norm_snw_vec = (df[ts, "SnowWater"] .- mean(df[ts, "SnowWater"])) ./ std(df[ts, "SnowWater"])
norm_slw_vec = (df[ts, "SoilWater"] .- mean(df[ts, "SoilWater"])) ./ std(df[ts, "SoilWater"])
nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)
m50_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
```

## Step 4: Pretrain the NNs
```julia
ep_grad_opt = HydroModels.GradOptimizer(component=ep_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
ep_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(ep_nn_flux)]))')
ep_output = (log_evap_div_lday=log_evap_div_lday_vec,)

ep_opt_params, epnn_loss_df = ep_grad_opt(
    [ep_input_matrix], [ep_output],
    tunable_pas=ComponentVector(nn=(epnn=ep_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)
```

```text
Training... 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:10        
((nn = (epnn = Float32[-0.01974994, -0.64867365, 0.13971744, -0.038864814, 1.112919, 0.06545781, -0.83492047, -0.43565416, -0.0061696735, 0.13157892  …  0.29205686, 0.46278232, -3.2228746, 0.12062144, 0.06469498, -2.8806944, 0.2044914, -1.7033366, 0.2584572, -0.266739])), 1001×4 DataFrame
  Row │ iter   loss        time                     params
      │ Int64  Float64     DateTime                 Array…
──────┼───────────────────────────────────────────────────────────────────────────────
    1 │     1  9.49144     2024-11-21T10:10:57.861  Float32[0.229399, -0.413821, 0.3…
    2 │     2  8.53332     2024-11-21T10:10:57.867  Float32[0.239399, -0.423821, 0.3…
    3 │     3  7.77423     2024-11-21T10:10:57.873  Float32[0.249313, -0.43382, 0.33…
    4 │     4  7.52688     2024-11-21T10:10:57.878  Float32[0.259011, -0.443769, 0.3…
    5 │     5  7.49975     2024-11-21T10:10:57.883  Float32[0.267003, -0.452015, 0.3…
    6 │     6  7.48733     2024-11-21T10:10:57.888  Float32[0.27375, -0.459019, 0.30…
    7 │     7  7.47953     2024-11-21T10:10:57.894  Float32[0.279562, -0.465082, 0.3…
  ⋮   │   ⋮        ⋮                  ⋮                             ⋮
  995 │   995  0.00583641  2024-11-21T10:11:04.627  Float32[-0.0199612, -0.650014, 0…
  996 │   996  0.00582027  2024-11-21T10:11:04.631  Float32[-0.0199189, -0.649745, 0…
  997 │   997  0.00580418  2024-11-21T10:11:04.640  Float32[-0.0198766, -0.649477, 0…
  998 │   998  0.00578814  2024-11-21T10:11:04.645  Float32[-0.0198344, -0.649209, 0…
  999 │   999  0.00577214  2024-11-21T10:11:04.650  Float32[-0.0197921, -0.648941, 0…
 1000 │  1000  0.0057562   2024-11-21T10:11:04.657  Float32[-0.0197499, -0.648674, 0…
 1001 │  1000  0.0057562   2024-11-21T10:11:04.662  Float32[-0.0197499, -0.648674, 0…
                                                                      987 rows omitted)
```

```julia
q_grad_opt = HydroModels.GradOptimizer(component=q_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
q_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(q_nn_flux)]))')
q_output = (log_flow=log_flow_vec,)

q_opt_params, qnn_loss_df = q_grad_opt(
    [q_input_matrix], [q_output],
    tunable_pas=ComponentVector(nn=(qnn=q_nn_params,)),
    const_pas=ComponentVector(),
    return_loss_df=true
)
```

```text
Training... 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:12        
((nn = (qnn = Float32[0.42091706, -0.2777198, 0.30891904, 0.4412809, 0.5341712, -0.3962445, -0.575717, -0.6653889, -0.45068413, 0.7776514  …  0.1986277, 0.22517225, 0.16427116, -1.7506562, -0.16015986, -1.1870662, -1.7102151, 0.41655824, 0.35249898, 0.10037236])), 1001×4 DataFrame
  Row │ iter   loss        time                     params
      │ Int64  Float64     DateTime                 Array…
──────┼───────────────────────────────────────────────────────────────────────────────
    1 │     1  1.25913     2024-11-21T10:12:32.773  Float32[0.235686, -0.42516, 0.36…
    2 │     2  1.21188     2024-11-21T10:12:32.778  Float32[0.245686, -0.41516, 0.35…
    3 │     3  1.18657     2024-11-21T10:12:32.815  Float32[0.253736, -0.406389, 0.3…
    4 │     4  1.12638     2024-11-21T10:12:32.820  Float32[0.262394, -0.397315, 0.3…
    5 │     5  1.02448     2024-11-21T10:12:32.826  Float32[0.2712, -0.38864, 0.3437…
    6 │     6  0.903213    2024-11-21T10:12:32.831  Float32[0.279957, -0.379719, 0.3…
    7 │     7  0.77082     2024-11-21T10:12:32.837  Float32[0.289023, -0.370426, 0.3…
  ⋮   │   ⋮        ⋮                  ⋮                             ⋮
  995 │   995  0.00126523  2024-11-21T10:12:39.001  Float32[0.420713, -0.277343, 0.3…
  996 │   996  0.00125197  2024-11-21T10:12:39.009  Float32[0.420763, -0.27742, 0.30…
  997 │   997  0.00123855  2024-11-21T10:12:39.015  Float32[0.420802, -0.277493, 0.3…
  998 │   998  0.00122546  2024-11-21T10:12:39.020  Float32[0.42084, -0.277568, 0.30…
  999 │   999  0.00121264  2024-11-21T10:12:39.025  Float32[0.420878, -0.277644, 0.3…
 1000 │  1000  0.00120045  2024-11-21T10:12:39.030  Float32[0.420917, -0.27772, 0.30…
 1001 │  1000  0.00120045  2024-11-21T10:12:39.037  Float32[0.420917, -0.27772, 0.30…
                                                                      987 rows omitted)
```

## Step 5: Retrain the NNs in the Exp-Hydro model
```julia
m50_opt = HydroModels.GradOptimizer(component=m50_model, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=100)
config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
norm_pas = ComponentVector(
    snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
    snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
)
m50_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...)
)
m50_tunable_pas = ComponentVector(
    nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
)
m50_opt_params, m50_loss_df = m50_opt(
    [m50_input], [q_output],
    tunable_pas=m50_tunable_pas,
    const_pas=m50_const_pas,
    config=[config],
    return_loss_df=true
)
```

## Step 6: Plot the results
```julia
```
