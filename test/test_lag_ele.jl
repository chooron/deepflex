include("../src/HydroModels.jl")
using ModelingToolkit: @variables, @parameters
using ComponentArrays

@variables slowflow fastflow slowflow_lag fastflow_lag
@parameters x1 x2

lag_funcs = [
    LumpedHydro.LagFlux(slowflow => slowflow_lag, x1, LumpedHydro.uh_3_half),
    LumpedHydro.LagFlux(fastflow => fastflow_lag, x2, LumpedHydro.uh_3_half),
]

lag_ele = LumpedHydro.LagElement(:lag, lfuncs=lag_funcs);
input_matrix = [1 2 3 2 3 4 1 2;] # 1 2 3 2 3 4 1 2

input_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([[1 2 3 2 3 4 1 2; 1 2 3 2 3 4 1 2]], 5))
input_arr = permutedims(input_arr, (1, 3, 2))
# input_arr = input_arr[1:1, :, :]
#* input_num * ts len * node num
cv = ComponentVector(NamedTuple{Tuple([Symbol(:node_, i) for i in 1:5])}(repeat([(x1=2.39, x2=5.39)], 5)))

output = LumpedHydro.run_multi_fluxes(lag_ele, input=input_arr, params=cv)
# output[1,1,:]