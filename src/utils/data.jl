##* used for flux calculate

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

function extract_params(params::NamedTuple, ::Vector{Symbol})
    params
end

function extract_params(::Vector, param_names::Vector{Symbol})
    @assert length(param_names) == 0
    NamedTuple()
end

function process_output(output::Union{T,Vector{T}}, output_names::Symbol) where {T<:Number}
    # namedtuple([output_names], [max.(T(0.0), output)])
    namedtuple([output_names], [output])
end

function process_output(output::Union{Vector{T},Vector{Vector{T}}}, output_names::Vector{Symbol}) where {T<:Number}
    # namedtuple(output_names, [max.(T(0.0), o) for o in output])
    namedtuple(output_names, output)
end

function process_output(output::Matrix{T}, output_names::Symbol) where {T<:Number}
    tmp_output = vec(output)
    tmp_output = length(tmp_output) > 1 ? [tmp_output] : tmp_output
    # namedtuple([output_names], [max.(T(0.0), output)])
    namedtuple([output_names], tmp_output)
end

function process_output(output::Matrix{T}, output_names::Vector{Symbol}) where {T<:Number}
    # namedtuple(output_names, [max.(T(0.0), o) for o in output])
    output_list = []
    for idx in length(output_names)
        tmp_output = output[idx, :]
        tmp_output = length(tmp_output) > 1 ? tmp_output : first(tmp_output)
        push!(output_list, tmp_output)
    end
    namedtuple(output_names, output_list)
end