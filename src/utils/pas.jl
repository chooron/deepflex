"""
获取demo输入数据
"""
get_demo_pas(flux::AbstractFlux) = ComponentVector(
    params=NamedTuple{Tuple(flux.infos[:param])}(ones(length(flux.infos[:param]))),
    initstates=NamedTuple{Tuple(flux.infos[:state])}(ones(length(flux.infos[:state]))),
)

get_demo_pas(flux::AbstractNeuralFlux) = ComponentVector(
    nn=NamedTuple{Tuple(flux.infos[:nn])}(zeros(flux.nnios[:paramlen]))
)

get_demo_pas(bucket::AbstractBucket) = ComponentVector(
    params=NamedTuple{Tuple(bucket.infos[:param])}(ones(length(flux.infos[:param]))),
    initstates=NamedTuple{Tuple(bucket.infos[:state])}(ones(length(bucket.infos[:state]))),
    nn=NamedTuple{Tuple(bucket.infos[:nn])}([zeros(flux.nnios[:paramlen]) for flux in bucket.fluxes if flux isa AbstractNeuralFlux]),
)

get_demo_pas(route::AbstractRoute) = ComponentVector(
    params=NamedTuple{Tuple(route.infos[:param])}(ones(length(route.infos[:param]))),
    initstates=NamedTuple{Tuple(route.infos[:state])}(ones(length(route.infos[:state]))),
    nn=NamedTuple{Tuple(route.infos[:nn])}([zeros(flux.nnios[:paramlen]) for flux in route.fluxes if flux isa AbstractNeuralFlux]),
)