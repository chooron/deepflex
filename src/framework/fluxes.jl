"""
水文非状态计算公式，对应superflexpy的ParameterizeElement
"""
@kwdef mutable struct SimpleFlux <: AbstractFlux
    input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}}
    output_names::Vector{Symbol}
    parameters::Union{NamedTuple,Nothing} = nothing
    func::Function
end

mutable struct LuxNNFlux <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    func::Function
    parameters::Any
end

## build flux
## ----------------------------------------------------------------------
function SimpleFlux(
    input_names::Union{Vector{Symbol},Vector{Dict{Symbol,Symbol}}},
    output_names::Vector{Symbol};
    parameters::ComponentVector{T},
    func::Function
) where {T<:Number}
    if length(parameters) == 0
        parameters = nothing
    else
        parameters = (; parameters...)
    end
    SimpleFlux(input_names, output_names, parameters, func)
end

function LuxNNFlux(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    lux_model::Lux.AbstractExplicitLayer,
    seed::Integer=42
)
    rng = MersenneTwister()
    Random.seed!(rng, seed)
    ps, st = Lux.setup(rng, lux_model)
    func = (x, p) -> lux_model(x, p, st)[1]
    LuxNNFlux(input_names, output_names, func, ComponentArray(ps=ps, st=st))
end

## ----------------------------------------------------------------------
## callable function
function (flux::SimpleFlux)(input::ComponentVector{T}) where {T<:Number}
    tmp_input = (; input[flux.input_names]...)
    tmp_output = flux.func(tmp_input, flux.parameters)
    ComponentVector(; Dict(flux.output_names[i] => tmp_output[i] for i in 1:length(flux.output_names))...)
end

function (flux::LuxNNFlux)(input::ComponentVector{T}) where {T<:Number}
    x = hcat([input[nm] for nm in flux.input_names]...)'
    y_pred = flux.func(x, flux.parameters[:ps])
    if size(y_pred, 2) > 1
        output = ComponentVector(; Dict(k => y_pred[i, :] for (i, k) in enumerate(flux.output_names))...)
    else
        output = ComponentVector(; Dict(k => first(y_pred[i, :]) for (i, k) in enumerate(flux.output_names))...)
    end
    return output
end

## ----------------------------------------------------------------------

## *namedtuple type generation function for SimpleFlux
## ---------------------------------------------------------------------- 
function gen_namedtuple_type(input_names::Vector{Symbol}, dtype::Union{Type,TypeVar})
    Union{NamedTuple{tuple(input_names...),NTuple{length(input_names),dtype}},
        NamedTuple{tuple(input_names...),NTuple{length(input_names),Vector{dtype}}}}
end
## ----------------------------------------------------------------------

## *parameters getter setter
## ----------------------------------------------------------------------
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
## ----------------------------------------------------------------------

## *训练后，更新模型内部参数
function update!(flux::LuxNNFlux, tstate)
    flux.parameters = tstate.parameters
end
## ----------------------------------------------------------------------
