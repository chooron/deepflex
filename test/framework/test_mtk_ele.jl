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

    DeepFlex.MTKElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

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
sf_ele = Surface_ExpHydro(name=:SF);
se_ele = Soil_ExpHydro(name=:SE);
@btime model = DeepFlex.MTKUnit(name=:model, elements=[sf_ele, se_ele])

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = (snowwater=0.0, soilwater=1303.004248)

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:10000
input = (time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
@btime sol = DeepFlex.solve_prob(
    model,
    input=input,
    params=params,
    init_states=init_states
)

sol_u = hcat(sol.u...)
sol
states = namedtuple(collect(keys(init_states)), [sol_u[1, :], sol_u[2, :]])

@btime result = model(input, states, params)
# todo 根据想获取的数据进行临时计算