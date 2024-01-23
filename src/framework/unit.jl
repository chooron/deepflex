# @kwdef mutable struct Unit{T} <: Component where {T<:Number}
#     id::String

#     # model structure
#     stucture::AbstractGraph

#     # inner variables
#     fluxes::Dict{Symbol,Vector{T}} = Dict()

#     # attribute
#     input_names::Vector{Symbol} = []
# end
function build_unit(; unit::Unit, paraminfos::Vector{ParamInfo})
    parameters = Dict(p.name => p.value for p in paraminfos)

end

function set_parameters!(unit::Unit; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for idx in topological_sort(unit.structure)
        tmp_ele = get_prop(unit.structure, idx, :ele)
        if isa(tmp_ele, ParameterizedElement) | isa(tmp_ele, StateParameterizedElement)
            set_parameters!(tmp_ele, paraminfos=paraminfos)
        end
        if isa(tmp_ele, StateElement) | isa(tmp_ele, StateParameterizedElement)
            set_states!(tmp_ele, paraminfos=paraminfos)
        end
    end
end

function get_output(unit::Unit; input::Dict{Symbol,Vector{T}}) where {T<:Number}
    # initialize unit fluxes
    unit.fluxes = input
    # traversal of the directed graph
    for idx in topological_sort(unit.structure)
        tmp_ele = get_prop(unit.structure, idx, :ele)
        tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
        merge!(unit.fluxes, tmp_fluxes)
    end
    return unit.fluxes
end

