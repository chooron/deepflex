# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

model = DeepFlex.HydroNode(
    :exphydro_node,
    units=[DeepFlex.ExpHydro_Unit(name=:exphydro1), DeepFlex.ExpHydro_Unit(name=:exphydro2)],
    routes=namedtuple([:exphydro1, :exphydro2],
        [DeepFlex.ExpHydro_RouteElement(name=:exphydro1), DeepFlex.ExpHydro_RouteElement(name=:exphydro2)])
)
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:100

unit_params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_input = (time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
unit_init_states = (snowwater=0.0, soilwater=1303.004248)

input = (exphydro1_unit=unit_input, exphydro2_unit=unit_input)
params = (
    exphydro1_unit=unit_params, exphydro2_unit=unit_params,
    exphydro1_route=NamedTuple(), exphydro2_route=NamedTuple(),
    weights=(exphydro2=0.5, exphydro1=0.5),
)
init_states = (exphydro1_unit=unit_init_states, exphydro2_unit=unit_init_states)
results = model(input, params, init_states, step=false) # 分开计算和同时计算消耗的资源差不多