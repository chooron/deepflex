@reexport module M100

using ..HydroModels
using ..HydroModels.Lux
using ..HydroModels: step_func

function M100_ELE(; name::Symbol, mtk::Bool=true)
    ann = Lux.Chain(
        Lux.Dense(4 => 16, Lux.tanh),
        Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 5, Lux.leakyrelu)
    )

    funcs = [
        NeuralFlux([:snowwater, :soilwater, :temp, :prcp] => [:log_et_lday, :log_q, :asinh_melt, :asinh_ps, :asinh_pr], :nn => ann),
        SimpleFlux([:snowwater, :asinh_melt] => [:melt], flux_funcs=[(i, p) -> @.[relu(step_func(i[1]) * sinh(i[2]))]]),
        SimpleFlux([:snowwater, :temp, :melt, :asinh_ps] => [:new_snowwater],
            flux_funcs=[(i, p) -> @.[i[1] + relu(sinh(i[4]) * step_func(-i[2]) - i[3])]]
        ),
        SimpleFlux([:soilwater, :lday, :melt, :log_et_lday, :log_q, :asinh_pr] => [:new_soilwater],
            flux_funcs=[(i, p) -> @.[i[1] + relu(sinh(i[6])) + i[3] - step_func(i[1]) * i[2] * exp(i[4]) - step_fct(i[1]) * exp(i[5])]],
        ),
    ]

    dfuncs = [
        StateFlux([:new_snowwater] => [:snowwater], funcs=funcs),
        StateFlux([:new_soilwater] => [:soilwater], funcs=funcs),
    ]

    HydroBucket(
        Symbol(name, :_ele_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end


function Model(; name::Symbol)
    build_unit(
        name=name,
        components=[M100_ELE(name=name, mtk=mtk)]
    )
end

end