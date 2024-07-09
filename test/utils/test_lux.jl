using Lux
using LuxCore
using StableRNGs
using Zygote

struct HydroLayer{F1} <: Lux.AbstractExplicitLayer
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    param_names::Vector{Symbol}
    layer_func::Function
    init_params::F1

    function HydroLayer(flux_names::Pair{Vector{Symbol},Vector{Symbol}}, param_names::Vector{Symbol}; layer_func::Function, init_params=Lux.zeros32)
        return new{typeof(init_params)}(
            flux_names[1],
            flux_names[2],
            param_names,
            layer_func,
            init_params
        )
    end
end

function Lux.initialparameters(rng, l::HydroLayer)
    return NamedTuple{Tuple(l.param_names)}(l.init_params(rng, len(l.param_names)))
end

Lux.initialstates(_, ::HydroLayer) = NamedTuple()

Lux.parameterlength(l::HydroLayer) = length(l.param_names)

Lux.statelength(::HydroLayer) = 0

function (l::HydroLayer)(x::AbstractMatrix, ps, st::NamedTuple)
    y = l.layer_func.(eachrow(x), Ref(ps))
    return y, st
end

function (l::HydroLayer)(x::AbstractVector, ps, st::NamedTuple)
    y = l.layer_func(x, ps)
    return y, st
end

tmp_layer = HydroLayer