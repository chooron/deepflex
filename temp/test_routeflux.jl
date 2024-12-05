using DataInterpolations
using OrdinaryDiffEq
using Zygote
using Plots

input_arr = [3 4 4 5 6 10 23 24 38 40 59 63 49 32 22 12 9 4]
input_itp = LinearInterpolation(input_arr[1, :], collect(1:length(input_arr)))



function discharge2!(du, u, p, t)
    s_river, q_in = u[1], u[2]
    q_out = (s_river+q_in) / (p[1] + 1)
    du[1] = q_in - q_out
    du[2] = q_out + input_itp(t) - q_in
end

prob = ODEProblem(discharge2!, [0.1, 3], (1, length(input_arr)), [0.3])
sol3 = solve(prob, Tsit5(), saveat=1.0)

function discharge2(u, p, t)
    s_river, q_in = u[1], u[2]
    @info [s_river q_in]
    q_out = (s_river+q_in) / (p[1] + 1)
    [s_river + q_in - q_out, q_out + input_itp(t)]
end

prob2 = DiscreteProblem(discharge2, [0.1, 3], (1, length(input_arr)), [0.3])
sol4 = solve(prob2, FunctionMap())

input_arr = [3 4 4 5 6 10 23 24 38 40 59 63 49 32 22 12 9 4]
s_river = zeros(length(input_arr))
q_in = zeros(length(input_arr))
q_out = zeros(length(input_arr))
s_river0, q_in0 = 0.1, 3
p = [0.3]
for i in eachindex(input_arr)
    q_out[i] = (s_river0) / (p[1] + 1)
    s_river[i] = s_river0 + q_in0 - q_out[i]
    q_in[i] = q_out[i] + input_arr[i]
    s_river0, q_in0 = s_river[i], q_in[i]
    @info [q_out[i] q_in[i] s_river[i]]
end


plot(sol3)
plot!(sol4)
savefig("test_routeflux.png")
