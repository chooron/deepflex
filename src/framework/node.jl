mutable struct Node{U<:Unit,T<:Number} <: Component
    id::String
    units::Dict{Symbol,U}
    weights::Dict{Symbol,T}
    area::T
    routing_func::Function{Tuple{Dict{Symbol,Vector{T}},Vector{Symbol}},Dict{Symbol,Vector{T}}}
    target_names::Vector{Symbol}
end

function Node(; id::String, units::Dict{Symbol,U}, weights::Dict{Symbol,T}, area::T, target_names::Vector{Symbol}) where {U<:Unit,T<:Number}
    Node{U,T}(id, units, weights, area, target_names)
end

function get_output(node::Node, input::Dict{Symbol,Vector{T}}) where {T<:Number}
    total_output::Dict{Symbol,Vector{T}} = Dict()
    # Sum the results of each unit with their respective weights
    for k in keys(node.units)
        tmp_ouput = get_output(node.units[k], input=input)
        for name in node.target_names
            if !haskey(total_output, name)
                total_output[name] = tmp_ouput[name] .* node.weights[k]
            else
                total_output[name] += tmp_ouput[name] .* node.weights[k]
            end
        end
    end
    total_output
end