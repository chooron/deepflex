mutable struct Node{T} <: Component where {T<:Number}
    id::String
    units::Dict{Symbol,Unit}
    weights::Dict{Symbol,T}
    area::T
    target_names::Vector{Symbol}
end

function get_output(node::Node, input::Dict{Symbol,Vector{T}}) where {T<:Number}
    total_output::Dict{Symbol, Vector{T}}
    for (n, w) in zip(node.units, node.weights)
        tmp_ouput = get_output(n, input=input)
        for name in node.target_names
            if haskey(total_output, name)
                total_output[name] = tmp_ouput[name].* w
            else
                total_output[name] += tmp_ouput[name].* w
            end
        end
    end
    total_output
end