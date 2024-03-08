function Summation(input_names::Vector{Symbol},
    output_names::Symbol;
    parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        output_names,
        parameters=parameters,
        func=(i::NamedTuple, p::NamedTuple, sf::Function) -> begin
            sum([i[k] for k in keys(i)])
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
        func=(i::NamedTuple, p::NamedTuple, sf::Function) -> begin
            tmp_input = i[first(keys(i))]
            [p[k] .* tmp_input for k in keys(p)]
        end
    )
end

"""
将多个通量无修改输出，可用于改名
"""
function Tranparent(input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol}=Symbol[])

    if length(output_names) == 0
        if input_names isa Vector{Symbol}
            output_names = input_names
        elseif input_names isa Vector{Dict{Symbol,Symbol}}
            output_names = collect(values(input_names))
        else
            @error "$(typeof(input_names)) is not available!"
        end
    end

    SimpleFlux(
        input_names,
        output_names,
        parameters=ComponentVector(),
        func=(i::NamedTuple, p::NamedTuple, sf::Function) -> [i[k] for k in keys(i)]
    )
end

function Differ(
    input_names::Dict{Symbol,Vector{Symbol}},
    output_names::Vector{Symbol};
)
    tmp_input_names = collect(union(map(x -> x, [Set(nms) for nms in values(input_names)])...))
    SimpleFlux(
        tmp_input_names,
        output_names,
        parameters=ComponentVector(),
        func=(i::NamedTuple, p::NamedTuple, sf::Function) -> [sum([i[nm] for nm in input_names[:In]]) - sum([i[nm] for nm in input_names[:Out]])]
    )
end