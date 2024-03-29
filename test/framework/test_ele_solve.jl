# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

function Surface_ExpHydro(; name::Symbol)
    funcs = [
        DeepFlex.PetFlux([:temp, :lday]),
        DeepFlex.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        DeepFlex.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DeepFlex.DifferFlux(Dict(:In => [:snowfall], :Out => [:melt]), :snowwater),
    ]

    DeepFlex.HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

ele = Surface_ExpHydro(name=:sf)

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = (snowwater=0.0, soilwater=1303.004248)

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:10000
input = (time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

@btime sol = DeepFlex.solve_prob(ele, input=input, params=params, init_states=init_states) # 1.21s
@btime sol = DeepFlex.solve_probv2(ele, input=input, params=params, init_states=init_states) # 36.292ms