#* used for flux calculate
function extract_input(input::Union{ComponentVector,NamedTuple}, input_names::Symbol)
    namedtuple([input_names], [input[input_names]])
end

function extract_input(input::Union{ComponentVector,NamedTuple}, input_names::Vector{Symbol})
    namedtuple(input_names, [input[k] for k in input_names])
end

function extract_input(input::Union{ComponentVector,NamedTuple}, input_names::Vector{Pair})
    namedtuple([k for (_, k) in input_names], [input[k] for (k, _) in input_names])
end

function extract_params(params::ComponentVector, param_names::Vector{Symbol})
    namedtuple(param_names, [params[k] for k in param_names])
end

function extract_params(params::NamedTuple, param_names::Vector{Symbol})
    params
end

function process_output(output_names::Symbol, output::Union{T,Vector{T}}) where {T<:Number}
    # namedtuple([output_names], [max.(T(0.0), output)])
    namedtuple([output_names], [output])
end

function process_output(output_names::Vector{Symbol}, output::Union{Vector{T},Vector{Vector{T}}}) where {T<:Number}
    # namedtuple(output_names, [max.(T(0.0), o) for o in output])
    ComponentVector(namedtuple(output_names, output))
end

function process_output(output_names::Symbol, output::Matrix{T}) where {T<:Number}
    # namedtuple([output_names], [max.(T(0.0), output)])
    ComponentVector(namedtuple([output_names], [vec(output)]))
end

function process_output(output_names::Vector{Symbol}, output::Matrix{T}) where {T<:Number}
    # namedtuple(output_names, [max.(T(0.0), o) for o in output])
    namedtuple(output_names, [output[idx, :] for idx in length(output_names)])
end