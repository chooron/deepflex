using Lux
using LuxCore
using Zygote
using StableRNGs
using ComponentArrays
using Optimization
using OptimizationOptimisers

# function LSTMCompact(in_dims, hidden_dims, out_dims)
#     lstm_cell = LSTMCell(in_dims => hidden_dims)
#     classifier = Dense(hidden_dims => out_dims, sigmoid)
#     return @compact(; lstm_cell, classifier) do x::AbstractArray{T,2} where {T}
#         x = reshape(x, size(x)..., 1)
#         x_init, x_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
#         y, carry = lstm_cell(x_init)
#         output = [vec(classifier(y))]
#         for x in x_rest
#             y, carry = lstm_cell((x, carry))
#             output = vcat(output, [vec(classifier(y))])
#         end
#         @return hcat(output...)
#     end
# end

function LSTMCompact(in_dims, hidden_dims, out_dims)
    lstm_cell = LSTMCell(in_dims => hidden_dims)
    classifier = Dense(hidden_dims => out_dims, sigmoid)
    return @compact(; lstm_cell, classifier) do x::AbstractArray{T, 2} where {T}
        x = reshape(x, size(x)..., 1)
        x_init, x_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
        y, carry = lstm_cell(x_init)
        output = [vec(classifier(y))]
        for x in x_rest
            y, carry = lstm_cell((x, carry))
            output = vcat(output, [vec(classifier(y))])
        end
        @return reduce(hcat, output)  # <--- This is the different line
    end
end

model = LSTMCompact(3, 10, 1)
ps, st = Lux.setup(StableRNGs.LehmerRNG(1234), model)
smodel = StatefulLuxLayer{true}(model, nothing, st)

# smodel(rand(3, 10), ComponentVector(ps))
# LuxCore.apply(smodel, rand(3, 10), ComponentVector(ps))
ps_axes = getaxes(ComponentVector(ps))
x = rand(3, 10)
y = rand(1, 10)
function object(u, p)
    ps = ComponentVector(u, ps_axes)
    sum((smodel(x, ps) .- y) .^ 2)
end

function objectv2(u, p)
    st = p[1]
    ps = ComponentVector(u, ps_axes)
    sum((model(x, ps, st)[1] .- y) .^ 2)
end

opt_func = Optimization.OptimizationFunction(objectv2, Optimization.AutoZygote())
opt_prob = Optimization.OptimizationProblem(opt_func, Vector(ComponentVector(ps)), [st])
opt_sol = Optimization.solve(opt_prob, OptimizationOptimisers.Adam(0.1), maxiters=100)


function MLPCompact(in_dims, hidden_dims, out_dims)
    mlp_cell = Dense(in_dims => hidden_dims)
    classifier = Dense(hidden_dims => out_dims, sigmoid)
    return @compact(; mlp_cell, classifier) do x::AbstractArray{T,2} where {T}
        output = classifier(mlp_cell(x))
        @return output
    end
end
