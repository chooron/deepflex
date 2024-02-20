mutable struct NNFlux <: AbstractFlux
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    func::Any
    parameters::Any
end

function NNFlux(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    model,
    seed=42
)
    rng = MersenneTwister()
    Random.seed!(rng, seed)
    ps, st = Lux.setup(rng, model)
    func = (x, p) -> model(x, p, st)[1]
    NNFlux(input_names, output_names, func, ComponentArray(ps=ps, st=st))
end


function update_lux_element!(ele::NNFlux, tstate)
    ele.model = tstate.model
    ele.parameters = tstate.parameters
    ele.states = tstate.states
end

function get_output(nn::NNFlux; input::ComponentVector{T}) where {T<:Number}
    x = hcat([input[nm] for nm in nn.input_names]...)'
    y_pred = nn.func(x, nn.parameters[:ps])
    if size(y_pred, 2) > 1
        output = ComponentVector(; Dict(k => y_pred[i, :] for (i, k) in enumerate(nn.output_names))...)
    else
        output = ComponentVector(; Dict(k => first(y_pred[i, :]) for (i, k) in enumerate(nn.output_names))...)
    end
    return output
end