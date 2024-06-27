function StdMeanNormFlux(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol}=[Symbol(:norm_, nm) for nm in input_names];
    param_names::Vector{Symbol}=vcat([[Symbol(:mean_, input_names), Symbol(:std_, input_names)]]...)
)
    inputs = [@variables $nm(t) = 0.0 for nm in input_names]
    outputs = [@variables $nm(t) = 0.0 for nm in output_names]
    params = [@parameters $nm = 0.0 for nm in param_names]
    params_ntp = NamedTuple{Tuple(param_names)}(params)
    flux_exprs = @.[(input - params_ntp[Symbol(:mean_, nm)]) / params_ntp[Symbol(:std_, nm)] for (input, nm) in zip(inputs, param_names)]

    SimpleFlux(input_names => output_names, params, flux_exprs=flux_exprs)
end

function MinMaxNormFlux(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol}=[Symbol(:norm_, nm) for nm in input_names],
    param_names::Vector{Symbol}=vcat([[Symbol(:max_, input_names), Symbol(:min_, input_names)]]...)
)
    SimpleFlux(input_names => output_names, param_names)
end

function StdMeanNormFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
    param_names::Vector{Symbol}=vcat([[Symbol(:mean_, input_names), Symbol(:std_, input_names)]]...)
)
    input_names, output_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm(t) = 0.0 for nm in input_names]
    outputs = [@variables $nm(t) = 0.0 for nm in output_names]
    params = [@parameters $nm = 0.0 for nm in param_names]
    params_ntp = NamedTuple{Tuple(param_names)}(params)
    flux_exprs = @.[(input - params_ntp[Symbol(:mean_, nm)]) / params_ntp[Symbol(:std_, nm)] for (input, nm) in zip(inputs, param_names)]

    SimpleFlux(inputs => outputs, params, flux_exprs=flux_exprs)
end

function MinMaxNormFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
    param_names::Vector{Symbol}=vcat([[Symbol(:max_, input_names), Symbol(:min_, input_names)]]...)
)
    input_names, output_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm(t) = 0.0 for nm in input_names]
    outputs = [@variables $nm(t) = 0.0 for nm in output_names]
    params = [@parameters $nm = 0.0 for nm in param_names]
    params_ntp = NamedTuple{Tuple(param_names)}(params)
    flux_exprs = @.[(input - params_ntp[Symbol(:min_, nm)]) / (params_ntp[Symbol(:max_, nm)] - params_ntp[Symbol(:min_, nm)])
               for (input, nm) in zip(inputs, param_names)]

    SimpleFlux(inputs => outputs, params, flux_exprs=flux_exprs)
end

function TranparentFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
)
    old_names, new_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm(t) = 0.0 for nm in old_names]
    outputs = [@variables $nm(t) = 0.0 for nm in new_names]
    flux_exprs = @.[var for var in inputs]

    SimpleFlux(inputs => outputs, Num[], flux_exprs=flux_exprs)
end