mutable struct LuxNN <: LuxNNFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    func::Any
    parameters::Any
end

function LuxNN(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    model,
    seed=42
)
    rng = MersenneTwister()
    Random.seed!(rng, seed)
    ps, st = Lux.setup(rng, model)
    func = (x, p) -> model(x, p, st)[1]
    LuxNN(input_names, output_names, func, ComponentArray(ps=ps, st=st))
end

function pretrain!(nn::LuxNNFunc; input::ComponentVector{T}, train_config...) where {T<:Number}
    x = hcat([input[nm] for nm in nn.input_names]...)
    y = hcat([input[nm] for nm in nn.output_names]...)'

    function prep_pred_NN_pretrain(model_, input_)
        (params) -> model_(input_, params)
    end

    pred_NN_pretrain_fct = prep_pred_NN_pretrain(nn.func, permutedims(x))

    function loss_NN_pretrain(params, batch)
        sum((pred_NN_pretrain_fct(params)' .- batch) .^ 2)
    end

    optf = Optimization.OptimizationFunction((θ, p) -> loss_NN_pretrain(θ, y), Optimization.AutoZygote())
    optprob = Optimization.OptimizationProblem(optf, nn.parameters[:ps])
    sol = Optimization.solve(optprob, Adam(0.01), maxiters=100)
    nn.parameters = ComponentArray(nn.parameters; Dict(:ps => sol.u)...)
end


function LinearNN(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    hidd_size::Int=32,
    hidd_layer::Int=1,
    seed=42
)
    input_layer = Dense(length(input_names), hidd_size)
    output_layer = Dense(hidd_size, length(output_names))
    hidd_layer = [Dense(hidd_size, hidd_size) for _ in 1:hidd_layer]
    model = Chain(input_layer, hidd_layer..., output_layer)
    LuxNN(input_names, output_names, model=model, seed=seed)
end

function update_lux_element!(ele::LuxNN, tstate)
    ele.model = tstate.model
    ele.parameters = tstate.parameters
    ele.states = tstate.states
end

function get_output(nn::LuxNN; input::ComponentVector{T}) where {T<:Number}
    x = hcat([input[nm] for nm in nn.input_names]...)'
    y_pred = nn.func(x, nn.parameters[:ps])
    if size(y_pred, 2) > 1
        output = ComponentVector(; Dict(k => y_pred[i, :] for (i, k) in enumerate(nn.output_names))...)
    else
        output = ComponentVector(; Dict(k => first(y_pred[i, :]) for (i, k) in enumerate(nn.output_names))...)
    end
    return output
end