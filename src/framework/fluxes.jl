
function get_func_infos(funcs::Vector{F}) where {F<:AbstractFlux}
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    parameter_names = Vector{Symbol}()

    for func in funcs
        if func.input_names isa Dict
            union!(input_names, Vector(keys(func.input_names)))
        else
            union!(input_names, func.input_names)
        end
        if func.output_names isa Vector
            union!(output_names, func.output_names)
        else
            union!(output_names, [func.output_names])
        end
        union!(parameter_names, func.parameter_names)
    end

    input_names, output_names, parameter_names
end

function get_d_func_infos(funcs::Vector{F}) where {F<:AbstractFlux}
    input_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    parameter_names = Vector{Symbol}()

    for func in funcs
        union!(state_names, func.output_names)
        union!(parameter_names, func.parameter_names)

        if func.input_names isa Dict
            for k in keys(func.input_names)
                union!(input_names, Vector(func.input_names[k]))
            end
        else
            union!(input_names, func.input_names)
        end
    end
    input_names, state_names, parameter_names
end

"""
Simple flux
"""
struct SimpleFlux <: AbstractFlux
    # attribute
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}}
    output_names::Union{Symbol,Vector{Symbol}}
    parameter_names::Vector{Symbol}
    # function 
    func::Function
    # step function
    step_func::Function
end

function SimpleFlux(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Union{Symbol,Vector{Symbol}},
    parameter_names::Vector{Symbol};
    func::Function
)
    SimpleFlux(
        input_names,
        output_names,
        parameter_names,
        func,
        DEFAULT_SMOOTHER
    )
end

mutable struct LuxNNFlux <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    func::Function
    parameters::Any
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
function (flux::SimpleFlux)(input::ComponentVector{T}, parameters::ComponentVector{T}) where {T<:Number}
    tmp_input = extract_input(input, flux.input_names)
    func_output = flux.func(tmp_input, NamedTuple(parameters[flux.parameter_names]), flux.step_func)
    process_output(flux.output_names, func_output)
end

function extract_input(input::ComponentVector{T}, input_names::Vector{Symbol}) where {T<:Number}
    namedtuple(input_names, [input[k] for k in input_names])
end

function extract_input(input::ComponentVector{T}, input_names::Dict{Symbol,Symbol}) where {T<:Number}
    namedtuple(collect(values(input_names)), [input[k] for k in keys(input_names)])
end

function process_output(output_names::Symbol, output::Union{T,Vector{T}}) where {T<:Number}
    ComponentVector(namedtuple([output_names], [output]))
end

function process_output(output_names::Vector{Symbol}, output::Union{Vector{T},Vector{Vector{T}}}) where {T<:Number}
    ComponentVector(namedtuple(output_names, output))
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
