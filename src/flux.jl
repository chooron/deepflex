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

function get_input_names(func::AF) where {AF<:AbstractFlux}
    if eltype(func.input_names) isa Pair
        input_names = [v for (_, v) in func.input_names]
    elseif func.input_names isa Vector
        input_names = func.input_names
    else
        input_names = [func.input_names]
    end
    input_names
end

function get_output_names(func::AF) where {AF<:AbstractFlux}
    if func.output_names isa Vector
        output_names = func.output_names
    else
        output_names = [func.output_names]
    end
    output_names
end

function get_param_names(func::AF) where {AF<:AbstractFlux}
    if func.param_names isa Vector
        param_names = func.param_names
    else
        param_names = [func.param_names]
    end
    param_names
end

function get_func_infos(funcs::Vector)
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    param_names = Vector{Symbol}()
    for func in funcs
        # extract the input and output name in flux
        tmp_input_names = get_input_names(func)
        tmp_output_names = get_output_names(func)
        tmp_param_names = get_param_names(func)
        # 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(tmp_output_names, input_names)
        # 输入需要排除已有的输出变量，表明这个输入是中间计算得到的
        tmp_input_names = setdiff(tmp_input_names, output_names)
        # 合并名称
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(param_names, tmp_param_names)
    end
    input_names, output_names, param_names
end

function get_dfunc_infos(dfuncs::Vector)
    input_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    param_names = Vector{Symbol}()

    for dfunc in dfuncs
        union!(input_names, get_input_names(dfunc))
        union!(state_names, get_output_names(dfunc))
        union!(param_names, dfunc.param_names)
    end
    input_names, state_names, param_names
end

## ----------------------------------------------------------------------
## callable function
function (flux::SimpleFlux)(input::NamedTuple, params::Union{ComponentVector,NamedTuple})
    tmp_input = extract_input(input, flux.input_names)
    tmp_params = extract_params(params, flux.param_names)
    func_output = flux.func(tmp_input, tmp_params, flux.step_func)
    process_output(flux.output_names, func_output)
end

function extract_input(input::NamedTuple, input_names::Symbol)
    namedtuple([input_names], [input[input_names]])
end

function extract_input(input::NamedTuple, input_names::Vector{Symbol})
    namedtuple(input_names, [input[k] for k in input_names])
end

function extract_input(input::NamedTuple, input_names::Vector{Pair})
    namedtuple([k for (_, k) in input_names], [input[k] for (k, _) in input_names])
end

function extract_params(params::ComponentVector, param_names::Vector{Symbol})
    namedtuple(param_names, [params[k] for k in param_names])
end

function extract_params(params::NamedTuple, param_names::Vector{Symbol})
    params
end

function process_output(output_names::Symbol, output::Union{T,Vector{T}}) where {T<:Number}
    namedtuple([output_names], [output])
end

function process_output(output_names::Vector{Symbol}, output::Union{Vector{T},Vector{Vector{T}}}) where {T<:Number}
    namedtuple(output_names, output)
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

function (flux::LagFlux)(input::NamedTuple, params::ComponentVector)
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
    sol = solve(prob)
    sol_u = hcat((sol.u)...)
    namedtuple([flux.input_names], [sol_u[1, :] .* l_input])
end

# TODO 关于Lux nn还需要进一步调整
struct NeuralFlux <: AbstractNNFlux
    input_names::Vector{Symbol}
    output_names::Union{Vector{Symbol},Symbol}
    param_names::Symbol
    func::Function
    nn::ODESystem
end

function NeuralFlux(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    param_names::Symbol,
    chain::Lux.AbstractExplicitLayer,
    seed::Int=42,
)
    func = (x, p) -> stateless_apply(model, x, p)[1]
    @named nn = NeuralNetworkBlock(length(input_names), length(output_names); chain=chain, rng=StableRNG(seed))
    NeuralFlux(input_names, output_names, param_names, func, nn)
end

function (flux::NeuralFlux)(input::NamedTuple, params::Union{ComponentVector,NamedTuple})
    x = hcat([input[nm] for nm in flux.input_names]...)'
    y_pred = flux.func(x, params[flux.param_names])
    if size(y_pred, 2) > 1
        output = namedtuple(flux.output_names, [y_pred[i, :] for (i, k) in enumerate(flux.output_names)])
    else
        output = namedtuple(flux.output_names, [first(y_pred[i, :]) for (i, k) in enumerate(flux.output_names)])
    end
    return output
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