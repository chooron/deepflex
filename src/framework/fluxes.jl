function get_input_names(func::AbstractFlux)
    if func.input_names isa Dict
        input_names = Vector(keys(func.input_names))
    else
        input_names = func.input_names
    end
    input_names
end

function get_output_names(func::AbstractFlux)
    if func.output_names isa Vector
        output_names = func.output_names
    else
        output_names = [func.output_names]
    end
    output_names
end

function get_func_infos(funcs::Vector{F}) where {F<:AbstractFlux}
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for func in funcs
        # extract the input and output name in flux
        tmp_input_names = get_input_names(func)
        tmp_output_names = get_output_names(func)
        # 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(tmp_output_names, tmp_input_names)
        # 输入需要排除已有的输出变量，表明这个输入是中间计算得到的
        tmp_input_names = setdiff(tmp_input_names, output_names)
        # 合并名称
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(param_names, func.param_names)
    end
    input_names, output_names, param_names
end

function get_dfunc_infos(funcs::Vector{F}) where {F<:AbstractFlux}
    input_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for func in funcs
        union!(input_names, get_input_names(func))
        union!(state_names, get_output_names(func))
        union!(param_names, func.param_names)
    end
    input_names, state_names, param_names
end

"""
Simple flux
"""
struct SimpleFlux <: AbstractFlux
    # attribute
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}}
    output_names::Union{Symbol,Vector{Symbol}}
    param_names::Vector{Symbol}
    # function 
    func::Function
    # step function
    step_func::Function
end

function SimpleFlux(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Union{Symbol,Vector{Symbol}};
    param_names::Vector{Symbol},
    func::Function
)
    SimpleFlux(
        input_names,
        output_names,
        param_names,
        func,
        DEFAULT_SMOOTHER
    )
end

struct LagFlux <: AbstractFlux
    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    param_names::Vector{Symbol}
    # function
    lag_times::Dict{Symbol,Symbol}
    lag_funcs::Dict{Symbol,Function}
    # step function
    step_func::Function

    # inner variable
    lag_states::ComponentVector
end

function LagFlux(
    lag_funcs::Dict{Symbol,Function};
    lag_times::Dict{Symbol,Symbol},
)
    input_names = collect(keys(lag_funcs))
    param_names = Vector{Symbol}()

    for value in lag_times
        union!(param_names, [value])
    end

    LagFlux(
        input_names,
        input_names,
        param_names,
        lag_times,
        lag_funcs,
        DEFAULT_SMOOTHER,
        ComponentVector()
    )
end

function init!(flux::LagFlux; parameters::ComponentVector{T}) where {T<:Number}
    # todo 添加dt插入功能
    lag_states = [init_lag_state(flux.lag_funcs[nm], parameters[flux.lag_times[nm]], 1.0) for nm in flux.input_names]
    ComponentVector(flux.input_names, lag_states)
end

function (flux::LagFlux)(input::ComponentVector{T}, parameters::ComponentVector{T}) where {T<:Number}
    tmp_output = T[]
    for nm in flux.input_names
        push!(tmp_output, lag_states[nm][1, 1] * input[nm] + lag_states[nm][2, 1])
        update_lag_state!(lag_states[nm], input[nm])
    end
    ComponentVector(namedtuple(flux.input_names, tmp_output))
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
    func_output = flux.func(tmp_input, NamedTuple(parameters[flux.param_names]), flux.step_func)
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

## ----------------------------------------------------------------------

## *训练后，更新模型内部参数
function update!(flux::LuxNNFlux, tstate)
    flux.parameters = tstate.parameters
end
## ----------------------------------------------------------------------
