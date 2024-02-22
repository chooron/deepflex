"""
水文非状态计算公式，对应superflexpy的ParameterizeElement
"""
@kwdef mutable struct SimpleFlux <: AbstractFlux
    input_names::Vector{Symbol}
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

mutable struct RoutingFlux <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    lag_states::ComponentVector
    lag_weights::ComponentVector
end

## build flux
## ----------------------------------------------------------------------
function SimpleFlux(input_names::Vector{Symbol},
    output_names::Vector{Symbol},
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

function RoutingFlux(
    input_names::Vector{Symbol};
    lag_time::Union{T,ComponentVector{T}},
    lag_func::Function
) where {T<:Number}

    if typeof(lag_time) == T
        lag_time = ComponentVector(; Dict(nm => lag_time for nm in input_names)...)
    end

    # init lag states
    lag_states = ComponentVector(; Dict(nm => begin
        zeros(Int(ceil(lag_time[nm])))
    end for nm in input_names)...)

    # build weight
    lag_weights = ComponentVector(; Dict(k => begin
        [
            lag_func(i + 1, ceil(lag_time[k])) - lag_func(i, ceil(lag_time[k]))
            for i in 1:(ceil(lag_time[k])|>Int)
        ]
    end for k in keys(lag_time))...)

    RoutingFlux(
        input_names,
        input_names,
        lag_states,
        lag_weights)
end
## ----------------------------------------------------------------------
## callable function
function (flux::SimpleFlux)(input::ComponentVector{T}) where {T<:Number}
    tmp_input = (; input[flux.input_names]...)
    flux.func(tmp_input, flux.parameters)
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


function (flux::RoutingFlux)(input::ComponentVector{T}) where {T<:Number}
    input = input[flux.input_names]
    solved_state = solve_lag(flux, input=input)

    # Get the new lag value to restart
    final_states = ComponentVector(; Dict(k => begin
        tmp_state = solved_state[k][end, :]
        tmp_state[:, 1:end-1] = tmp_state[:, 2:end]
        tmp_state[:, end] = T(0)
    end for k in keys(flux.input_names))...)

    flux.lag_states = final_states

    output = ComponentVector(; Dict(k => solved_state[k][:, 1] for k in flux.input_names)...)
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

function solve_lag(flux::RoutingFlux; input::ComponentVector{T}) where {T<:Number}
    max_weight_len = maximum([length(flux.lag_weights[k]) for k in keys(flux.lag_weights)])
    max_input_len = maximum([length(input[k]) for k in keys(input)])
    output = ComponentVector(; Dict(k => zeros(T, max_input_len, max_weight_len) for k in keys(input))...)

    for k in keys(output)
        for (w, ls, i) in zip(flux.lag_weights[k], flux.lag_states[k], input[k])
            for ts in 1:max_input_len
                updated_state = ls .+ i[ts] .* w
                output[k][ts, 1:length(w)] .= updated_state
                @info updated_state
                ls = vcat(updated_state[2:end], 0)
            end
        end
    end
    return output
end
