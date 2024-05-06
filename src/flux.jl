"""
Simple flux
"""
struct SimpleFlux <: AbstractFlux
    # attribute
    input_names::Union{Symbol,Vector{Symbol},Vector{Pair}}
    output_names::Union{Symbol,Vector{Symbol}}
    param_names::Vector{Symbol}
    # function 
    func::Function
    # step function
    step_func::Function

    function SimpleFlux(
        input_names::Union{Symbol,Vector{Symbol},Vector{Pair}},
        output_names::Union{Symbol,Vector{Symbol}};
        param_names::Vector{Symbol},
        func::Function
    )
        return new(
            input_names,
            output_names,
            param_names,
            func,
            DEFAULT_SMOOTHER
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
    func=(i::NamedTuple, p::NamedTuple, sf::Function) ->
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
    func=(i::NamedTuple, p::NamedTuple, sf::Function) ->
        @.((i[input_names] - p[param_names[2]]) / (p[param_names[1]] - p[param_names[2]]))
)

DifferFlux(
    influx_names::Vector{Symbol},
    outflux_names::Vector{Symbol},
    state_names::Symbol;
) = SimpleFlux(
    vcat(influx_names, outflux_names),
    state_names,
    param_names=Symbol[],
    func=(i::NamedTuple, p::NamedTuple, sf::Function) -> sum([i[nm] for nm in influx_names]) -
                                                         sum([i[nm] for nm in outflux_names])
)

function (flux::SimpleFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, flux.input_names)
    tmp_params = extract_params(params, flux.param_names)
    func_output = flux.func(tmp_input, tmp_params, flux.step_func)
    process_output(flux.output_names, func_output)
end

struct LagFlux <: AbstractFlux
    # * 这里设计成只针对一个变量的Flux
    # attribute
    input_names::Symbol
    output_names::Symbol
    param_names::Symbol
    # function
    lag_func::Function
    # step function
    step_func::Function
end

function LagFlux(
    input_names::Symbol,
    output_names::Symbol;
    lag_func::Function,
    param_names::Symbol,
    step_func::Function=DEFAULT_SMOOTHER
)
    LagFlux(
        input_names,
        output_names,
        param_names,
        lag_func,
        step_func,
    )
end

function (flux::LagFlux)(input::Union{ComponentVector,NamedTuple}, params::Union{ComponentVector,NamedTuple})
    l_input = input[flux.input_names]
    #* 首先将lagflux转换为discrete problem
    function discrete_prob!(du, u, p, t)
        tmp_u = l_input[Int(t)] .* p + u
        tmp_u = circshift(tmp_u, -1)
        tmp_u[end] = 0
        du[1:end] = tmp_u
    end
    lag_time = params[flux.param_names]
    delta_t = 1.0
    ts = 1:(ceil(lag_time / delta_t)|>Int)
    #* 将weight作为param输入到prob中
    lag_weights = [flux.lag_func(t, lag_time, flux.step_func) for t in ts]
    lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
    prob = DiscreteProblem(discrete_prob!, zeros(eltype(l_input), length(ts)), (1.0, length(l_input)), lag_weights)
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    sol_u = hcat((sol.u)...)
    namedtuple([flux.input_names], [sol_u[1, :] .* l_input])
end


struct NeuralFlux <: AbstractNNFlux
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
    process_output(flux.output_names, y_pred)
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