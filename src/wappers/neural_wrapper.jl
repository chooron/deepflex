struct NeuralWrapper{N,F,M} <: AbstractNeuralWrapper
    model::N
    func::F
    meta::M

    function NeuralWrapper(fluxes::Pair, model::N, model_name::Symbol; rng::R=StableRNG(42)) where {N,R}
        #* assert the chain has a name
        #* extract the parameter and state
        ps, st = Lux.setup(rng, model)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func(x, p) = LuxCore.apply(model, x, ComponentVector(p, ps_axes), st)[1]
        meta = HydroMeta(name=Symbol(model_name, :_nwrapper), inputs=fluxes.first, outputs=fluxes.second, nn_names=[model_name])
        return new{N,typeof(nn_func),typeof(meta)}(model, nn_func, meta)
    end
end

function (wrapper::NeuralWrapper)(input::Array{T,2}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nn][wrapper.meta.name]
    output = wrapper.func(input, nn_ps)
    return output
end

function (wrapper::NeuralWrapper)(input::Array{T,3}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nn][wrapper.meta.name]
    output = wrapper.func.(Matrix.(eachslice(input, dims=2)), Ref(nn_ps))
    return permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (1, 3, 2))
end