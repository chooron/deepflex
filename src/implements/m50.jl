@reexport module M50

using ..LumpedHydro
using ..LumpedHydro.NamedTupleTools
import ..LumpedHydro: Lux

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        LumpedHydro.PetFlux([:temp, :lday]),
        LumpedHydro.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        LumpedHydro.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        LumpedHydro.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        LumpedHydro.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        LumpedHydro.StateFlux([:snowfall], [:melt], :snowwater),
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Soil(; name::Symbol, mtk::Bool=true)

    # 神经网络的定义是在模型之内，需要提取到模型的参数
    et_ann = Lux.Chain(
        Lux.Dense(3 => 16, Lux.tanh),
        # Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu)
    )

    q_ann = Lux.Chain(
        Lux.Dense(2 => 16, Lux.tanh),
        # Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu)
    )

    funcs = [
        # normalize
        LumpedHydro.StdNormFlux(:snowwater, :norm_snw),
        LumpedHydro.StdNormFlux(:soilwater, :norm_slw),
        LumpedHydro.StdNormFlux(:temp, :norm_temp),
        LumpedHydro.StdNormFlux(:prcp, :norm_prcp),
        # ET ANN
        LumpedHydro.NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, param_names=:etnn, chain=et_ann),
        # Q ANN
        LumpedHydro.NeuralFlux([:norm_slw, :norm_prcp], :flow, param_names=:qnn, chain=q_ann),
        # 一些变量转为用于状态计算的形式
        LumpedHydro.SimpleFlux([:soilwater, :lday, :evap], :real_evap, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:soilwater] * i[:lday] * i[:evap])),
        LumpedHydro.SimpleFlux([:soilwater, :flow], :real_flow, param_names=Symbol[],
            func=(i, p; kw...) -> kw[:smoooth_func](i[:soilwater]) * exp(i[:flow])),
    ]

    dfuncs = [
        LumpedHydro.StateFlux([:infiltration], [:real_evap, :real_flow], :soilwater)
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroUnit(
        name,
        surface=Surface(name=name, mtk=mtk),
        soil=Soil(name=name, mtk=mtk),
        freewater=FreeWater(name=name),
        step=step,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        LumpedHydro.SimpleFlux(:flow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:flow])
    ]

    LumpedHydro.HydroElement(
        Symbol(name, :_route_),
        funcs=funcs,
        mtk=mtk
    )
end

"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""

function Node(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroNode(
        name,
        units=[Unit(name=name, mtk=mtk, step=step)],
        routes=[Route(name=name)],
    )
end


end