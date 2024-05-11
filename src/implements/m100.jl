@reexport module M100

function Model(; name::Symbol, mtk::Bool=true)
    ann = Lux.Chain(
        Lux.Dense(4 => 16, Lux.tanh),
        # Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 4, Lux.leakyrelu)
    )

    DeepFlex.NeuralFlux([:SnowWater, :SoilWater, :Temp, :Prcp],
        [:Snowfall, :Rainfall, :Melt, :Evap, :Flow],
        param_names=:etnn, chain=ann
    )

end


function Node(; name::Symbol)
    elements = [
    ]
    build_unit(name=name, elements=elements)
end

end