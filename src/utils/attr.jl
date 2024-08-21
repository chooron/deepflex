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

#* used for getting element attr
get_ode_func(::AbstractElement) = nothing
get_ode_func(buc::AbstractBucket) = buc.ode_func
