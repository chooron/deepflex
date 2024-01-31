mutable struct LuxNN <: AbstractFunc
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    model::Any
    device::Any
    parameters::Any
    states::Any
end

function LuxNN(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    model,
    device=cpu_device(),
    seed=42
)
    rng = MersenneTwister()
    Random.seed!(rng, seed)
    ps, st = Lux.setup(rng, model) .|> device
    LuxNN(input_names, output_names, model, device, ps, st)
end


function LinearNN(
    input_names::Vector{Symbol},
    output_names::Vector{Symbol};
    hidd_size::Int=32,
    hidd_layer::Int=1,
    device=cpu_device(),
    seed=42
)
    input_layer = Dense(length(input_names), hidd_size)
    output_layer = Dense(hidd_size, length(output_names))
    hidd_layer = [Dense(hidd_size, hidd_size) for _ in 1:hidd_layer]
    model = Chain(input_layer, hidd_layer..., output_layer)
    LuxNN(input_names, output_names, model=model, device=device, seed=seed)
end

function update_lux_element!(ele::LuxNN, tstate)
    ele.model = tstate.model
    ele.parameters = tstate.parameters
    ele.states = tstate.states
end

function get_output(ele::LuxNN; input::ComponentVector{T}) where {T<:Number}
    x = hcat([input[k] for k in keys(ele.input_names)]...)'
    y_pred = cpu_device()(Lux.apply(ele.model, ele.device(x), ele.parameters, ele.states)[1])
    y_pred = dropdims(y_pred, dims=1)
    output = ComponentVector(; Dict(k => y_pred[i] for (i, k) in enumerate(ele.output_names))...)
    return output
end