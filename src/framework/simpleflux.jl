"""
水文非状态计算公式，对应superflexpy的ParameterizeElement
"""
struct SimpleFlux <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    parameters::Union{NamedTuple,Nothing}
    func::Function
end

function build_flux(input_names::Vector{Symbol},
    output_names::Vector{Symbol},
    parameters::Union{ComponentVector{T},Nothing},
    func::Function) where {T<:Number}

    if parameters !== nothing
        parameters = (; parameters...)
    end

    SimpleFlux(
        input_names,
        output_names,
        parameters,
        func
    )
end

function (flux::SimpleFlux)(input::ComponentVector{T}) where {T<:Number}
    tmp_input = (; input[flux.input_names]...)
    flux.func(tmp_input, flux.parameters)
end

function gen_namedtuple_type(input_names::Vector{Symbol}, dtype::Union{Type,TypeVar})
    Union{NamedTuple{tuple(input_names...),NTuple{length(input_names),dtype}},
        NamedTuple{tuple(input_names...),NTuple{length(input_names),Vector{dtype}}}}
end

function get_parameters(func::SimpleFlux; names::Vector{Symbol}=nothing)::Dict{Symbol,Any}
    if isnothing(names)
        return func.parameters
    else
        return Dict(name => func.parameters[name] for name in names)
    end
end

function set_parameters!(func::SimpleFlux; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        if p.name in keys(func.parameters)
            func.parameters[p.name] = p.value
        end
    end
end