"""
Simple flux
"""
struct SimpleFlux <: AbstractSimpleFlux
    # attribute
    input_names::Union{Symbol,Vector{Symbol},Vector{Pair}}
    output_names::Union{Symbol,Vector{Symbol}}
    param_names::Vector{Symbol}
    # function 
    func::Function
    # step function
    smooth_func::Function

    function SimpleFlux(
        input_names::Union{Symbol,Vector{Symbol},Vector{Pair}},
        output_names::Union{Symbol,Vector{Symbol}};
        param_names::Vector{Symbol},
        func::Function,
        smooth_func::Function=step_func,
    )
        return new(
            input_names,
            output_names,
            param_names,
            func,
            smooth_func
        )
    end
end

StdNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:mean_, input_names), Symbol(:std_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple, kwargs...) ->
        @.((i[input_names] - p[param_names[1]]) / p[param_names[2]])
)

MinMaxNormFlux(
    input_names::Symbol,
    output_names::Symbol=Symbol(:norm_, input_names);
    param_names::Vector{Symbol}=[Symbol(:max_, input_names), Symbol(:min_, input_names)]
) = SimpleFlux(
    input_names,
    output_names,
    param_names=param_names,
    func=(i::NamedTuple, p::NamedTuple, kwargs...) ->
        @.((i[input_names] - p[param_names[2]]) / (p[param_names[1]] - p[param_names[2]]))
)

function (flux::AbstractSimpleFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    tmp_output = flux.func(tmp_input, tmp_params, smooth_func=flux.smooth_func)
    process_output(tmp_output, get_output_names(flux))
end

mutable struct StateFlux <: AbstractStateFlux
    influx_names::Vector{Symbol}
    outflux_names::Vector{Symbol}
    state_names::Symbol
    params_names::Vector{Symbol}
    func::Function

    function StateFlux(
        influx_names::Vector{Symbol},
        outflux_names::Vector{Symbol},
        state_names::Symbol;
        params_names::Vector{Symbol}=Symbol[],
        func::Function=(i::NamedTuple, p::NamedTuple; kwargs...) ->
            sum([i[nm] for nm in kwargs[:influx_names]]) - sum([i[nm] for nm in kwargs[:outflux_names]])
    )
        return new(
            influx_names,
            outflux_names,
            state_names,
            params_names,
            func
        )
    end
end

function (flux::AbstractStateFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    tmp_output = flux.func(tmp_input, tmp_params, influx_names=flux.influx_names, outflux_names=flux.outflux_names)
    process_output(tmp_output, get_output_names(flux))
end

struct LagFlux <: AbstractLagFlux
    # * 这里设计成只针对一个变量的Flux
    # attribute
    flux_name::Symbol
    # lag time name
    lag_time::Symbol
    # function
    lag_func::Function
    smooth_func::Function

    function LagFlux(
        flux_name::Symbol;
        lag_time::Symbol,
        lag_func::Function,
        smooth_func::Function=step_func,
    )
        new(
            flux_name,
            lag_time,
            lag_func,
            smooth_func
        )
    end
end

function solve_lag_weights(flux::LagFlux; input::NamedTuple, params::Union{ComponentVector,NamedTuple})
    #* 首先将lagflux转换为discrete problem
    function discrete_prob!(du, u, p, t)
        tmp_u = input[Int(t)] .* p + u
        tmp_u = circshift(tmp_u, -1)
        tmp_u[end] = 0
        du[1:end] = tmp_u
    end
    lag_time = params[flux.param_names]
    delta_t = 1.0
    ts = 1:(ceil(lag_time / delta_t)|>Int)
    #* 将weight作为param输入到prob中
    lag_weights = [flux.lag_func(t, lag_time, smooth_func=flux.smooth_func) for t in ts]
    lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
    prob = DiscreteProblem(discrete_prob!, zeros(eltype(input), length(ts)), (1.0, length(input)), lag_weights)
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    sol_u = hcat((sol.u)...)
    #* 得到权重计算结果
    sol_u[1, :]
end

function (flux::LagFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, get_input_names(flux))
    tmp_params = extract_params(params, get_param_names(flux))
    lag_weight = solve_lag_weights(flux, input=tmp_input, params=tmp_params)
    process_output(lag_weight .* tmp_input[get_input_names(flux)], get_output_names(flux))
end

struct NeuralFlux <: AbstractNeuralFlux
    input_names::Vector{Symbol}
    output_names::Union{Vector{Symbol},Symbol}
    param_names::Symbol
    func::Function
    nn::ODESystem
end

function NeuralFlux(
    input_names::Vector{Symbol},
    output_names::Union{Symbol,Vector{Symbol}};
    param_names::Symbol,
    chain::Lux.AbstractExplicitLayer,
    seed::Int=42,
)
    func = (x, p) -> LuxCore.stateless_apply(chain, x, p)
    nn = NeuralNetworkBlock(length(input_names), length(output_names); chain=chain, rng=StableRNG(seed), name=param_names)
    NeuralFlux(input_names, output_names, param_names, func, nn)
end

function (flux::NeuralFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    x = hcat([input[nm] for nm in flux.input_names]...)'
    y_pred = flux.func(x, params[flux.param_names])
    process_output(y_pred, get_output_names(flux))
end

# macro simpleflux(input_names, output_names, param_names)
#     func_nm = Symbol(output_names, :_func)
#     println(func_nm)
#     quote
#         simple_flux = SimpleFlux(
#             $(esc(input_names)),
#             $(esc(output_names)),
#             param_names=$(esc(param_names)),
#             func=$(func_nm),
#         )
#         simple_flux
#     end
# end