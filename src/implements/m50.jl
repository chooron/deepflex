@reexport module M50
using ..LumpedHydro
using ..NamedTupleTools
import ..LumpedHydro: Lux

"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""
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

    et_ann = Lux.Chain(
        Lux.Dense(3 => 16, Lux.tanh),
        Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu)
    )

    q_ann = Lux.Chain(
        Lux.Dense(2 => 16, Lux.tanh),
        Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu)
    )

    funcs = [
        # normalize
        StdMeanNormFlux([:snowwater, :soilwater, :prcp, :temp] => [:norm_snw, :norm_slw, :norm_prcp, :norm_temp],
            [:snowwater => [:mean_snowwater, :std_snowwater], :soilwater => [:mean_soilwater, :std_soilwater],
                :prcp => [:mean_prcp, :std_prcp], :temp => [:mean_temp, :std_temp]]),
        # ET ANN
        NeuralFlux([:norm_snw, :norm_slw, :norm_temp] => [:evap], :etnn => et_ann),
        # Q ANN
        NeuralFlux([:norm_slw, :norm_prcp] => [:flow], :qnn => q_ann),
        SimpleFlux([:soilwater, :lday, :evap] => [:realevap], flux_funcs=[(i, p) -> [i[1] * i[2] * i[3]]]),
        SimpleFlux([:soilwater, :flow] => [:realflow], flux_funcs=[(i, p) -> [ifelse(i[1] > 0, 1.0, 0.0) * exp(i[2])]]),
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