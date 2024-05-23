StdMeanNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:mean_, input_names), Symbol(:std_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple; kw...) ->
        @.[(i[input_names] - p[param_names[1]]) / p[param_names[2]]]
)

MinMaxNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:max_, input_names), Symbol(:min_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple; kw...) ->
        @.[(i[input_names] - p[param_names[2]]) / (p[param_names[1]] - p[param_names[2]])]
)

TranparentFlux(
    old_names::Symbol,
    new_names::Symbol,
) = SimpleFlux(
    old_names,
    new_names,
    param_names=Symbol[],
    func=(i::NamedTuple, p::NamedTuple; kw...) -> [i[old_names]]
)

TranparentFlux(
    old_names::Vector{Symbol},
    new_names::Vector{Symbol},
) = SimpleFlux(
    old_names,
    new_names,
    param_names=Symbol[],
    func=(i::NamedTuple, p::NamedTuple; kw...) -> [i[nm] for nm in old_names]
)