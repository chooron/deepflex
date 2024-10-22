using DataInterpolations
using OrdinaryDiffEq
using Zygote

input_arr = [3 4 4 5 6 10 23 24 38 40 59 63 49 32 22 12 9 4]
input_itp = LinearInterpolation(input_arr[1, :], collect(1:length(input_arr)))
function discharge!(du, u, p, t)
    s_river, q_in = u[1], u[2]
    q_out = (s_river + q_in) / (p[1] + 1)
    du[1] = q_in - q_out
end

function condition(u, t, integrator)
    return true
end

function affect!(integrator)
    q_in = integrator.u[2]
    q_out = (integrator.uprev[1] + q_in) / (p[1] + 1)
    integrator.u[2] = q_out + input_itp(integrator.t) - q_in
end

cb = DiscreteCallback(condition, affect!)
prob = ODEProblem(discharge!, [0.1, 3], (1, length(input_arr)), [p])
sol = solve(prob, Tsit5(), dt=1.0, callback=cb)
sol_u = Array(sol)


integrator = init(DiscreteProblem(discharge!, [0.1, 3], (1, length(input_arr)), [0.3]), FunctionMap(), saveat=1.0)
for _ in 1:length(input_arr)
    step!(integrator)
    q_out = (integrator.u[1] + integrator.u[2]) / (integrator.p[1] + 1)
    integrator.u[2] = q_out + input_itp(integrator.t-1) - integrator.u[2]
end

sol2 = integrator.sol
# Zygote.gradient(loss1, 0.2)


function discharge2!(du, u, p, t)
    s_river, q_in = u[1], u[2]
    q_out = (s_river + q_in) / (p[1] + 1)
    du[1] = q_in - q_out
    du[2] = q_out + input_itp(t) - q_in
end

prob = DiscreteProblem(discharge2!, [0.1, 3], (1, length(input_arr)), [0.3])
sol2 = solve(prob, FunctionMap(), saveat=1.0)


function discharge3!(du, u, p, t)
    s_river, q_in = u[1], u[2]
    q_out = (s_river + q_in) / (p[1] + 1)
    du[1] = q_in - q_out
    du[2] = q_out + input_itp(t) - q_in
end

s_river = zeros(length(input_arr))
q_in = zeros(length(input_arr))
q_out = zeros(length(input_arr))
s_river0, q_in0 = 0.1, 3
p = [0.3]
for i in eachindex(input_arr)
    q_out[i] = (s_river0 + q_in0) / (p[1] + 1)
    s_river[i] = q_in0 - q_out[i]
    q_in[i] = q_out[i] + input_arr[i] - q_in0
    @info [q_out[i] q_in[i] s_river[i]]
    s_river0, q_in0 = s_river[i], q_in[i]
end


