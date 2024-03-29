# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

function Surface_GR4J(; name::Symbol)
    funcs = [
        DeepFlex.RainfallFlux([:prcp, :pet]),
        DeepFlex.SimpleFlux([:prcp, :pet], :pet,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(sf(i[:pet] - i[:prcp]) * (i[:pet] - i[:prcp]))),
        DeepFlex.InfiltrationFlux([:rainfall])
    ]

    DeepFlex.MTKElement(
        name=name,
        funcs=funcs,
        dfuncs=[],
    )
end

function Soil_GR4J(; name::Symbol)

    funcs = [
        DeepFlex.SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        DeepFlex.EvapFlux([:soilwater, :pet], param_names=[:x1]),
        DeepFlex.PercolationFlux([:soilwater], param_names=[:x1]),
    ]

    dfuncs = [
        DeepFlex.Differ(Dict(:In => [:infiltration], :Out => [:evap, :percolation]), :soilwater)
    ]

    DeepFlex.MTKElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function Routing_GR4J(; name::Symbol)

    funcs = [
        DeepFlex.SimpleFlux([:infiltration, :percolation, :saturation], :tempflow, param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        DeepFlex.SimpleFlux([:tempflow], [:tempslowflow], param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.9)
        DeepFlex.SimpleFlux([:tempflow], [:tempfastflow], param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.1)
        DeepFlex.LagFlux([:tempslowflow], [:slowflow], param_names=[:x4], lag_func=uh_1_half),
        DeepFlex.LagFlux([:tempfastflow], [:fastflow], param_names=[:x4], lag_func=uh_2_full),
        DeepFlex.RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        DeepFlex.SimpleFlux([:routingstore], :routedflow, param_names=[:x3, :γ],
            func=(i, p, sf) -> @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:routingstore]^p[:γ]))),
        DeepFlex.Summation([:routedflow, :recharge, :fastflow], :flow)
    ]

    dfuncs = [
        DeepFlex.DifferFlux(Dict(:In => [:slowflow, :recharge], :Out => [:routedflow]), :routingstore),
    ]

    DeepFlex.MTKElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end