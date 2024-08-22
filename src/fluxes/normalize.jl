function MinMaxNormFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
    params::Vector{Vector{Num}}
)
    flux_exprs = [(input - param[2]) / (param[1] - param[2]) for (input, param) in zip(fluxes[1], params)]
    params_vec = reduce(vcat, params)

    SimpleFlux(
        fluxes,
        params_vec,
        exprs=flux_exprs,
    )
end

function StdMeanNormFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
    params::Vector{Vector{Num}}
)
    flux_exprs = [(input - param[1]) / param[2] for (input, param) in zip(fluxes[1], params)]
    params_vec = reduce(vcat, params)

    SimpleFlux(
        fluxes,
        params_vec,
        exprs=flux_exprs,
    )
end

function TranparentFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
)
    flux_exprs = [var for var in fluxes[1]]
    SimpleFlux(fluxes, Num[], exprs=flux_exprs)
end

function StdMeanNormFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
    param_names::Vector{Symbol}
)
    input_names, output_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm for nm in input_names]
    outputs = [@variables $nm for nm in output_names]
    params = [first(@parameters $names[1:2]) for names in param_names]

    flux_exprs = [(input - param[1]) / param[2] for (input, param) in zip(inputs, params)]

    function inner_flux_func(input::AbstractVector, params::AbstractVector)
        [(i - p[1]) / p[2] for (i, p) in zip(input, params)]
    end

    function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
        reduce(hcat, inner_flux_func.(eachrow(input), Ref(params)))
    end

    SimpleFlux(
        NamedTuple{(Tuple(input_names))}(inputs),
        NamedTuple{(Tuple(output_names))}(outputs),
        NamedTuple{(Tuple(param_names))}(params),
        inner_flux_func,
        flux_exprs,
    )
end

function MinMaxNormFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
    param_names::Vector{Symbol}
)
    input_names, output_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm = 0.0 for nm in input_names]
    outputs = [@variables $nm = 0.0 for nm in output_names]
    params = [first(@parameters $names[1:2]) for names in param_names]
    flux_exprs = [(input - param[2]) / (param[1] - param[2]) for (input, param) in zip(inputs, params)]

    function inner_flux_func(input::AbstractVector, params::AbstractVector)
        [(i - p[2]) / (p[1] - p[2]) for (i, p) in zip(input, params)]
    end

    function inner_flux_func(input::AbstractMatrix, params::AbstractVector)
        reduce(hcat, inner_flux_func.(eachrow(input), Ref(params)))
    end

    SimpleFlux(
        NamedTuple{(Tuple(input_names))}(inputs),
        NamedTuple{(Tuple(output_names))}(outputs),
        NamedTuple{(Tuple(param_names))}(params),
        inner_flux_func,
        flux_exprs,
    )
end

function TranparentFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
)
    old_names, new_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm = 0.0 for nm in old_names]
    outputs = [@variables $nm = 0.0 for nm in new_names]
    flux_exprs = [var for var in inputs]

    SimpleFlux(inputs => outputs, Num[], flux_exprs=flux_exprs)
end
