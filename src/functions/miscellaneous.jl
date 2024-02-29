function Summation(input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=(i, p) -> begin
            [sum([i[k] for k in keys(i)])]
        end
    )
end



"""
将一个通量按比例拆分
"""
function Splitter(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol};
    parameters::ComponentVector{T}) where {T<:Number}
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=(i::NamedTuple, p::NamedTuple) -> begin
            tmp_input = i[first(keys(i))]
            [p[k] .* tmp_input for k in keys(p)]
        end
    )
end

"""
将多个通量无修改输出，可用于改名
"""
function Tranparent(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=Symbol[];
    parameters::ComponentVector=ComponentVector())

    if length(output_names) == 0
        output_names = input_names
    end

    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=(i::NamedTuple, p::Nothing=nothing) -> [i[k] for k in keys(i)]
    )
end