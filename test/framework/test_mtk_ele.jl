# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
include("../../src/DeepFlex.jl")

function Soil_ExpHydro(; name::Symbol)
    funcs = [
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        DeepFlex.BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        DeepFlex.SurfaceflowFlux([:soilwater], param_names=[:Smax]),
        DeepFlex.FlowFlux([:baseflow, :surfaceflow])
    ]

    dfuncs = [
        DeepFlex.DifferFlux(Dict(:In => [:infiltration], :Out => [:evap, :flow]), :soilwater)
    ]

    DeepFlex.MTKElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

# build 
ele = Soil_ExpHydro(name=:se);

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowwater=0.0, soilwater=1303.004248)

file_path = "data/cache/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
input = ComponentVector(time=1:1000, pet=df[1:1000, :Pet], infiltration=df[1:1000, :Infiltration])

sol = DeepFlex.solve_prob(
    ele,
    input=input,
    params=params,
    init_states=init_states
)