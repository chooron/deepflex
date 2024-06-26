#! 测试 build_state_func
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using Symbolics

include("../../src/LumpedHydro.jl")

@variables routingstore(t) recharge(t) slowflow_lag(t) fastflow_lag(t) routedflow(t)
@parameters x2, x3, ω, γ
fluxes = [
    LumpedHydro.SimpleFlux([:routingstore] => [:recharge], [:x2, :x3, :ω]),
    LumpedHydro.SimpleFlux([:routingstore, :recharge, :slowflow_lag] => [:routedflow], [:x3, :γ]),
    LumpedHydro.SimpleFlux([:routedflow, :recharge, :fastflow_lag] => [:flow])
]

state_expr = sum([slowflow_lag]) - sum([routedflow, recharge])
state_expr = LumpedHydro.build_state_func(fluxes, state_expr, [:slowflow_lag, :routedflow, :recharge])