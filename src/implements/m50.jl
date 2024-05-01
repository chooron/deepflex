@reexport module M50

using ..DeepFlex
import ..DeepFlex: Lux

"""
SoilWaterReservoir in Exp-Hydro
"""
function Surface(; name::Symbol, mtk::Bool=true)
    funcs = [
        DeepFlex.PetFlux([:temp, :lday]),
        DeepFlex.SnowfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.MeltFlux([:snowwater, :temp], param_names=[:Tmax, :Df]),
        DeepFlex.RainfallFlux([:prcp, :temp], param_names=[:Tmin]),
        DeepFlex.InfiltrationFlux([:rainfall, :melt])
    ]

    dfuncs = [
        DeepFlex.DifferFlux([:snowfall], [:melt], :snowwater),
    ]

    DeepFlex.HydroElement(
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
        DeepFlex.StdNormFlux(:snowwater, :norm_snw),
        DeepFlex.StdNormFlux(:soilwater, :norm_slw),
        DeepFlex.StdNormFlux(:temp, :norm_temp),
        DeepFlex.StdNormFlux(:prcp, :norm_prcp),
        # ET ANN
        DeepFlex.NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, param_names=:etnn, chain=et_ann),
        # Q ANN
        DeepFlex.NeuralFlux([:norm_slw, :norm_prcp], :flow, param_names=:qnn, chain=q_ann),
    ]

    dfuncs = [
        DeepFlex.SimpleFlux([:soilwater, :infiltration, :lday, :evap, :flow], :soilwater,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:infiltration] -
                                  sf(i[:soilwater]) * i[:lday] * exp(i[:evap]) -
                                  sf(i[:soilwater]) * exp(i[:flow]))
        )
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        mtk=mtk
    )
end

function Route(; name::Symbol, mtk::Bool=true)

    funcs = [
        DeepFlex.SimpleFlux([:flow], :flow, param_names=Symbol[], func=(i, p, sf) -> i[:flow])
    ]

    DeepFlex.HydroElement(
        Symbol(name, :_route_),
        funcs=funcs,
        mtk=mtk
    )
end

"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""

function Node(; name::Symbol, mtk::Bool=true)
    elements = [
        Surface(name=name,mtk=mtk),
        Soil(name=name,mtk=mtk)
    ]
    DeepFlex.HydroNode(
        name,
        units=elements
    )
end




# function M100(; name::Symbol, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
#     elements = [
#         LinearNN(
#             [:SnowWater, :SoilWater, :Temp, :Prcp],
#             [:Snowfall, :Rainfall, :Melt, :Evap, :Flow],
#             hnamed_size=32,
#             hnamed_layer=1
#         ),
#         SnowWater_ExpHydro_ODE(
#             name=:sr,
#             parameters=parameters,
#             init_states=init_states[[:SnowWater]]
#         ),
#         SoilWater_M50_ODE(
#             name=:wr,
#             parameters=parameters,
#             init_states=init_states[[:SoilWater]]
#         )
#     ]
#     build_unit(name=name, elements=elements)
# end

end