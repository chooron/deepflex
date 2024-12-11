struct NamedTupleIOAdapter{C,M} <: AbstractIOAdapter where {C<:AbstractComponent,M<:HydroMeta}
    component::C
    meta::M

    function NamedTupleIOAdapter(component::C) where {C<:AbstractComponent}
        new{C,typeof(component.meta)}(component, component.meta)
    end
end

function (adapter::NamedTupleIOAdapter)(
    input::NamedTuple,
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
)::NamedTuple
    @assert((all(input_name in keys(input) for input_name in get_input_names(adapter))),
        "Missing required inputs. Expected all of $(get_input_names(adapter)), but got $(keys(input)).")
    input_matrix = Matrix(reduce(hcat, [input[k] for k in get_input_names(adapter)])')
    output_matrix = adapter.component(input_matrix, pas; config=config, kwargs...)
    return NamedTuple{Tuple(vcat(get_state_names(adapter), get_output_names(adapter)))}(eachslice(output_matrix, dims=1))
end

function (adapter::NamedTupleIOAdapter)(
    input::Vector{<:NamedTuple},
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...
)::Vector{<:NamedTuple}
    for i in eachindex(input)
        @assert(all(input_name in keys(input[i]) for input_name in get_input_names(adapter)),
            "Missing required inputs. Expected all of $(get_input_names(adapter)), but got $(keys(input[i])) at $i input.")
    end
    input_mats = [reduce(hcat, collect(input[i][k] for k in get_input_names(adapter))) for i in eachindex(input)]
    input_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), input_mats)
    output_arr = adapter.component(input_arr, pas; config=config, kwargs...)
    output_names_tuple = Tuple(vcat(get_state_names(adapter), get_output_names(adapter)))
    return [NamedTuple{output_names_tuple}(eachslice(output_arr_, dims=1)) for output_arr_ in eachslice(output_arr, dims=2)]
end