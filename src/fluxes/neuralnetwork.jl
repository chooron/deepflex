# todo  构建一些常用的神经网络flux

function MLPNNFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
    hidden_size::Integer,
    hidden_layers::Integer,
    activation::Function,
)
    Lux.Chain(
        Lux.Dense(hidden_size, hidden_size, activation),
        (i for i in 2:hidden_layers)...,
        Lux.Dense(hidden_size, 1),
    )

end