function Summation(input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters,
        (i, p) -> begin
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
        parameters,
        (i::NamedTuple, p::NamedTuple) -> begin
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
        parameters,
        (i::NamedTuple, p::NamedTuple) -> [i[k] for k in keys(i)]
    )
end

function Differ(
    input_names::Vector{Dict{Symbol,Vector{Symbol}}},
    output_names::Vector{Symbol};
)
    tmp_input_names = vcat([collect(nms) for nms in values(input_names)])
    SimpleFlux(
        tmp_input_names,
        output_names,
        parameters,
        (i::NamedTuple, p::NamedTuple) -> sum([i[nm] for nm in input_names[:In]]) - sum([i[nm] for nm in input_names[:Out]])
    )
end