@reexport module M50

"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""
using ..LumpedHydro
using ..NamedTupleTools
import ..LumpedHydro: Lux

function SurfaceStorage(; name::Symbol, mtk::Bool=true)
    funcs = [
        SimpleFlux([:temp, :lday] => [:pet]),
        SimpleFlux([:prcp, :temp] => [:snowfall], [:Tmin]),
        SimpleFlux([:snowwater, :temp] => [:melt], [:Tmax, :Df]),
        SimpleFlux([:prcp, :temp] => [:rainfall], [:Tmin]),
        SimpleFlux([:rainfall, :melt] => [:infiltration]),
    ]

    dfuncs = [
        StateFlux([:snowfall] => [:melt], :snowwater),
    ]

    HydroElement(
        Symbol(name, :_surf_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function SoilStorage(; name::Symbol, mtk::Bool=true)

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
        StdMeanNormFlux([:snowwater, :soilwater, :prcp, :temp] => [:norm_snw, :norm_slw, :norm_prcp, :norm_temp]),
        # ET ANN
        NeuralFlux([:norm_snw, :norm_slw, :norm_temp] => [:evap], :etnn => et_ann),
        # Q ANN
        NeuralFlux([:norm_slw, :norm_prcp] => [:flow], :qnn => q_ann),
        # 一些变量转为用于状态计算的形式
        SimpleFlux([:soilwater, :lday, :evap] => [:realevap], Symbol[], flux_funcs=[(i, p) -> @.(i[1] * i[2] * i[3])]),
        SimpleFlux([:soilwater, :flow] => [:realflow], Symbol[], flux_funcs=[(i, p) -> @.(ifelse_func(i[1]) * exp(i[2]))]),
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
        elements=[SurfaceStorage(name=name, mtk=mtk), SoilStorage(name=name, mtk=mtk)],
        step=step,
    )
end

end