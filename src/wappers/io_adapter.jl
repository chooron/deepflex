struct NamedTupleIOAdapter{N,M} <: AbstractIOAdapter where {N,M<:ComponentVector}
    component::AbstractComponent
    meta::M

    function NamedTupleIOAdapter(component::C; name::Union{Symbol,Nothing}=nothing) where {C<:AbstractComponent}
        @assert !(typeof(component) isa AbstractIOAdapter) "Component $component is already an IO adapter"
        wrapper_name = isnothing(name) ? Symbol("##wrapper#", get_name(component)) : name
        new{wrapper_name,typeof(component.meta)}(component, component.meta)
    end
end

function (adapter::NamedTupleIOAdapter)(
    input::NamedTuple,
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
)::NamedTuple
    for input_name in get_input_names(adapter.component)
        @assert(input_name in keys(input), "Missing required inputs. Expected $input_name in input, but got $(keys(input)).")
    end
    input_matrix = Matrix(reduce(hcat, [input[k] for k in get_input_names(adapter.component)])')
    output_matrix = adapter.component(input_matrix, pas; config=config, kwargs...)
    output_names_tuple = Tuple(vcat(get_state_names(adapter.component), get_output_names(adapter.component)))
    return NamedTuple{output_names_tuple}(eachslice(output_matrix, dims=1))
end

function (adapter::NamedTupleIOAdapter)(
    input::Vector{<:NamedTuple},
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
)::Vector{<:NamedTuple}
    for i in eachindex(input)
        for input_name in get_input_names(adapter.component)
            @assert(input_name in keys(input[i]), "Missing required inputs. Expected $input_name in input, but got $(keys(input)).")
        end
    end
    input_mats = [reduce(hcat, collect(input[i][k] for k in get_input_names(adapter.component))) for i in eachindex(input)]
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), input_mats), (2, 3, 1))
    output_arr = adapter.component(input_arr, pas; config=config, kwargs...)
    output_names_tuple = Tuple(vcat(get_state_names(adapter.component), get_output_names(adapter.component)))
    return [NamedTuple{output_names_tuple}(eachslice(output_arr_, dims=1)) for output_arr_ in eachslice(output_arr, dims=2)]
end