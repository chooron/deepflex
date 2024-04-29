using ModelingToolkitNeuralNets
using ModelingToolkit
import ModelingToolkit.t_nounits as t
import ModelingToolkit.D_nounits as Dt
using ModelingToolkitStandardLibrary.Blocks
using OrdinaryDiffEq
using Optimization
using OptimizationOptimisers: Adam
using SciMLStructures
using SciMLStructures: Tunable
using SymbolicIndexingInterface
using StableRNGs
using Lux
using Plots

Fbrk = 100.0
vbrk = 10.0
Fc = 80.0
vst = vbrk / 10
vcol = vbrk * sqrt(2)
function friction(v)
    sqrt(2 * MathConstants.e) * (Fbrk - Fc) * exp(-(v / vst)^2) * (v / vst) +
    Fc * tanh(v / vcol)
end

function friction_true()
    @variables y(t) = 0.0
    @constants Fu = 120.0
    eqs = [
        Dt(y) ~ Fu - friction(y)
    ]
    return ODESystem(eqs, t, name=:friction_true)
end

model_true = structural_simplify(friction_true())
prob_true = ODEProblem(model_true, [], (0, 0.1), [])
sol_ref = solve(prob_true, Rodas4(); saveat=0.001)

scatter(sol_ref, label="velocity")
scatter(sol_ref.t, friction.(first.(sol_ref.u)), label="friction force")

function friction_ude(Fu)
    @variables y(t) = 0.0
    @constants Fu = Fu
    @named nn_in = RealInputArray(nin=2)
    @named nn_out = RealOutputArray(nout=1)
    eqs = [Dt(y) ~ Fu - nn_out.u[1]
        y ~ nn_in.u[1]
        y ~ nn_in.u[2]
        ]
    return ODESystem(eqs, t, name=:friction, systems=[nn_in, nn_out]) # 
end

Fu = 120.0
model = friction_ude(Fu)

chain = Lux.Chain(
    Lux.Dense(2 => 10, Lux.tanh, use_bias=false),
    Lux.Dense(10 => 10, Lux.mish, use_bias=false),
    Lux.Dense(10 => 1, use_bias=false)
)
@named nn = NeuralNetworkBlock(2, 1; chain=chain, rng=StableRNG(1111))

# eqs = [connect(model.nn_in, nn.output)
#     connect(model.nn_out, nn.input)]

eqs = [connect(model.nn_in, nn.input)
    connect(model.nn_out, nn.output)]

ude_sys = complete(ODESystem(eqs, t, systems=[model, nn], name=:ude_sys))
sys2 = structural_simplify(ude_sys)
for p in ModelingToolkit.defaults(sys2)
    # if typeof(p[2]) <: AbstractVector
    println(p)
    # end
end
prob = ODEProblem(sys2, [], (0, 1), Float32[])
solve(prob, Rodas4())
# get_vars = getu(sys, [sys.friction.y])
# get_refs = getu(model_true, [model_true.y])
# x0 = reduce(vcat, getindex.((default_values(sys),), tunable_parameters(sys)))

# function loss(x, (prob, sol_ref, get_vars, get_refs))
#     new_p = SciMLStructures.replace(Tunable(), prob.p, x)
#     new_prob = remake(prob, p = new_p, u0 = eltype(x).(prob.u0))
#     ts = sol_ref.t
#     new_sol = solve(new_prob, Rodas4(), saveat = ts, abstol = 1e-8, reltol = 1e-8)
#     loss = zero(eltype(x))
#     for i in eachindex(new_sol.u)
#         loss += sum(abs2.(get_vars(new_sol, i) .- get_refs(sol_ref, i)))
#     end
#     if SciMLBase.successful_retcode(new_sol)
#         loss
#     else
#         Inf
#     end
# end

# cb = (opt_state, loss) -> begin
#     @info "step $(opt_state.iter), loss: $loss"
#     return false
# end

# of = OptimizationFunction{true}(loss, AutoForwardDiff())
# op = OptimizationProblem(of, x0, (prob, sol_ref, get_vars, get_refs))
# res = solve(op, Adam(5e-3); maxiters = 1000, callback = cb)

# res_p = SciMLStructures.replace(Tunable(), prob.p, res)
# res_prob = remake(prob, p = res_p)
# res_sol = solve(res_prob, Rodas4(), saveat = sol_ref.t)

# initial_sol = solve(prob, Rodas4(), saveat = sol_ref.t)

# scatter(sol_ref, idxs = [model_true.y], label = "ground truth velocity")
# plot!(res_sol, idxs = [sys.friction.y], label = "velocity after training")
# plot!(initial_sol, idxs = [sys.friction.y], label = "velocity before training")

# scatter(sol_ref.t, friction.(first.(sol_ref.u)), label = "ground truth friction")
# plot!(res_sol.t, Fu .- first.(res_sol(res_sol.t, Val{1}).u),
#     label = "friction from neural network")