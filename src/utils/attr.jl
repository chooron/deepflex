#* get variables for element
get_input_vars(flux::AbstractFlux) = flux.inputs

get_output_vars(flux::AbstractFlux) = flux.outputs
get_output_vars(::AbstractStateFlux) = Num[]

get_state_vars(::AbstractFlux) = Num[]
get_state_vars(flux::AbstractStateFlux) = [flux.state]

get_param_vars(flux::AbstractFlux) = flux.params
get_param_vars(::AbstractNeuralFlux) = Num[]

get_nnparam_vars(::AbstractFlux) = Symbolics.Arr[]
get_nnparam_vars(flux::AbstractNeuralFlux) = [flux.nnparam]

get_all_vars(flux::AbstractFlux) = reduce(union, [get_input_vars(flux), get_output_vars(flux), get_state_vars(flux)])

get_exprs(flux::AbstractFlux) = flux.exprs
get_exprs(flux::AbstractStateFlux) = [flux.expr]
get_exprs(flux::AbstractNeuralFlux) = [flux.expr]

#* used for getting element attr
get_ode_func(::AbstractElement) = nothing
get_ode_func(buc::AbstractBucket) = buc.ode_func

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