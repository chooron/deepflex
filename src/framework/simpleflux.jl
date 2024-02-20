"""
水文非状态计算公式，对应superflexpy的ParameterizeElement
"""
struct SimpleFlux{T<:Number} <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    parameters::Union{ComponentVector{T},Nothing}
    func::Function
end

function SimpleFlux(input_names::Vector{Symbol}, output_names::Vector{Symbol}, parameters::Union{ComponentVector{T},Nothing}, func::Function) where {T<:Number}
    SimpleFlux{T}(
        input_names,
        output_names,
        parameters,
        func
    )
end

function (flux::SimpleFlux)(input::ComponentVector{T}) where {T<:Number}
    tmp_input = input[flux.input_names]
    flux.func(tmp_input, flux.parameters)
end

function get_output(flux::SimpleFlux, input::ComponentVector{T}) where {T<:Number}
    output = Dict{Symbol,Vector{T}}(k => [] for k in flux.output_names)
    for idx in 1:length(input[first(keys(input))])
        tmp_input = ComponentVector(; Dict(k => input[k][idx] for k in flux.input_names)...)
        tmp_output = flux.func(tmp_input, flux.parameters)
        for k in flux.output_names
            push!(output[k], tmp_output[k])
        end
    end
    ComponentVector(; output...)
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