#* used for flux calculate
function extract_input(input::NamedTuple, input_names::Symbol)
    namedtuple([input_names], [input[input_names]])
end

function extract_input(input::NamedTuple, input_names::Vector{Symbol})
    namedtuple(input_names, [input[k] for k in input_names])
end

function extract_input(input::NamedTuple, input_names::Vector{Pair})
    namedtuple([k for (_, k) in input_names], [input[k] for (k, _) in input_names])
end

function extract_params(params::ComponentVector, param_names::Vector{Symbol})
    namedtuple(param_names, [params[k] for k in param_names])
end

function extract_params(params::NamedTuple, param_names::Vector{Symbol})
    params
end

function process_output(output_names::Symbol, output::Union{T,Vector{T}}) where {T<:Number}
    namedtuple([output_names], [output])
end

function process_output(output_names::Vector{Symbol}, output::Union{Vector{T},Vector{Vector{T}}}) where {T<:Number}
    namedtuple(output_names, output)
end