@reexport module M50

using ..LumpedHydro
using ..NamedTupleTools
import ..LumpedHydro: Lux

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        PetFlux([:temp, :lday]),
        SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        StateFlux([:snowfall], [:melt], :snowwater),
    ]

    HydroElement(
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
        StdMeanNormFlux(:snowwater, :norm_snw),
        StdMeanNormFlux(:soilwater, :norm_slw),
        StdMeanNormFlux(:temp, :norm_temp),
        StdMeanNormFlux(:prcp, :norm_prcp),
        # ET ANN
        NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, chain_name=:etnn, chain=et_ann),
        # Q ANN
        NeuralFlux([:norm_slw, :norm_prcp], :totalflow, chain_name=:qnn, chain=q_ann),
        # 一些变量转为用于状态计算的形式
        SimpleFlux([:soilwater, :lday, :evap], :realevap, param_names=Symbol[],
            func=(i, p; kw...) -> @.(i[:soilwater] * i[:lday] * i[:evap])),
        SimpleFlux([:soilwater, :totalflow], :realflow, param_names=Symbol[],
            func=(i, p; kw...) -> begin
                sf = kw[:smooth_func]
                @.(sf(i[:soilwater]) * exp(i[:totalflow]))
            end),
    ]

    dfuncs = [
        StateFlux([:infiltration], [:realevap, :realflow], :soilwater)
    ]

    HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Unit(; name::Symbol, mtk::Bool=true, step::Bool=true)
    HydroUnit(
        name,
        elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk)],
        step=step,
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        SimpleFlux(:totalflow, :flow, param_names=Symbol[], func=(i, p; kw...) -> i[:totalflow])
    ]

    HydroElement(
        name,
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