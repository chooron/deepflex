"""
Regional PRNN model

Jiang S., Zheng Y., & Solomatine D.. (2020) Improving AI system awareness of geoscience knowledge:
Symbiotic integration of physical approaches and deep learning. *Geophysical Research Letters*, 47. DOI: 10.1029/2020GL088229
"""
function RegionalPRNN(; name::String, parameters::ComponentVector{T}, init_states::ComponentVector{T}, solver::AbstractSolver=nothing) where {T<:Number}
    # attribute names
    attr_names = []

    hidd_size = parameters[:hidd_size]

    param_est_model = Lux.Chain(
        Lux.Dense(input_dim, hidd_size),
        Lux.Dense(hidd_size, 1, leakyrelu)
    )

    # paramter estimator
    region_param_est = NNEstimator(attr_names, [:f, :Smax, :Qmax, :Df, :Tmax, :Tmin], model=param_est_model)

    # 使用默认参数预先构建好ele
    # PRNN
    elements = [
        Surface_ExpHydro(
            name=:sf,
            init_states=init_states[[:SnowWater]]
        ),
        Soil_ExpHydro(
            name=:sl,
            init_states=init_states[[:SoilWater]]
        ),
        Routing_ExpHydro(name=:rt),
    ]

    input_dim = 4 + length(attr_names) # predict + forcing + attributes

    conv_model = Lux.Chain(
        Lux.Dense(input_dim, hidd_size, tanh),
        Lux.Dense(hidd_size, hidd_size, leakyrelu),
        Lux.Dense(hidd_size, 1, leakyrelu)
    )

    nn_input_nms = []
    funcs = [LuxNNFlux(nn_input_nms, [:Flow], lux_model=conv_model, seed=42)]
    push!(elements, SimpleElement(name=name, parameters=parameters, funcs=funcs))

    # 当存在estimator使需要重设参数
    build_unit(name=name, elements=elements, solver=solver, estimator=estimator)
end