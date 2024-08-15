#* hbv nn 耦合模块
module HBV_NN
#* 模型设计
#* 动态参数: 1. W, 使用(soiltype, landuse, soilwater)来预测 W
#* 动态参数: 2. L k0 k1 k2 kp, 使用(s1, s2, soiltype, landuse, prcpeff) 来预测
#* 后处理器: 3. 使用后处理神经网络通过q0,q1,q2来预测最终的输出结果
using ..LumpedHydro
using ..LumpedHydro: step_func
using Lux
using Symbolics: @variables
using ModelingToolkit: @parameters

"""
SoilWaterReservoir in HBV_EDU
"""
function SoilStorage(; name::Symbol)
    @variables soilwater prcpeff prcp pet evap W
    @parameters PWP

    funcs = [
        #* 直接用prcp替换原模型中的liquidwater
        NeuralFlux([soilwater, prcp] => [W], Lux.Chain(Lux.Dense(2 => 16), Lux.Dense(16 => 1, Lux.sigmoid_fast), name=:wnn)),
        SimpleFlux([soilwater, prcp, W] => [prcpeff], exprs=@.[prcp * W]),
        SimpleFlux([soilwater, pet] => [evap], [PWP], exprs=@.[(pet * min(1.0, soilwater / PWP))]),
    ]

    dfuncs = [
        StateFlux([prcp] => [prcpeff, evap], soilwater),
    ]

    HydroBucket(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function FreeWaterStorage(; name::Symbol)
    @variables s1 s2 q0 q1 q2 qp prcpeff flow k0 k1 soiltype landuse
    @parameters L k2 kp

    funcs = [
        NeuralFlux([s1, prcpeff] => [k0, k1], Lux.Chain(Lux.Dense(2 => 16), Lux.Dense(16 => 2, Lux.sigmoid_fast), name=:ksnn)),
        SimpleFlux([s1, k0] => [q0], [L], exprs=@.[max(0.0, s1 - L) * (k0 * 0.8 + 0.2)]),
        SimpleFlux([s1, k1] => [q1], exprs=@.[s1 * (k1 * 0.2 + 0.01)]),
        SimpleFlux([s1] => [qp], [kp], exprs=@.[s1 * kp]),
        SimpleFlux([s2] => [q2], [k2], exprs=@.[s2 * k2]),
        # todo 这个可以再根据q0, q1, q2构建预测模型
        # SimpleFlux([q0, q1, q2] => [flow], exprs=@.[q0 + q1 + q2]),
        NeuralFlux([s1, prcpeff, q0, q1, q2] => [flow], Lux.Chain(Lux.Dense(5 => 32), Lux.Dense(32 => 1, Lux.leakyrelu), name=:qnn)),
    ]

    dfuncs = [
        StateFlux([prcpeff] => [q0, q1, qp], s1),
        StateFlux([qp] => [q2], s2),
    ]

    HydroBucket(
        Symbol(name, :_zone_),
        funcs=funcs,
        dfuncs=dfuncs,
    )
end

function Unit(; name::Symbol)
    elements = [
        SoilStorage(name=name),
        FreeWaterStorage(name=name),
    ]

    HydroUnit(
        name,
        components=elements
    )
end
end