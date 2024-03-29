# *现在有个方式是将lagflux的计算过程转换成一个ode element
# *首先我们先用常规方式定义这种ode Element
# *这种element的state是通过lag time设置来定义的
using OrdinaryDiffEq
using DataInterpolations

x4 = 3.5
delta_t = 1.0
input = Float64[1, 2, 3, 4, 5, 3, 2, 1, 4]
itp = LinearInterpolation(input, 1:length(input))

function ode_sys!(du, u, p, t)
    tmp_u = itp(t) .* p + u
    tmp_u = circshift(tmp_u, -1)
    tmp_u[end] = 0
    println("t: $t, curr u: $tmp_u")
    du[1:end] = tmp_u
    # du[end] = t
    # print("$t" * " : ")
end

# !这个差值是始终在一个基值上的，即以1为基值，当值到2时，中间逐步累加的过程最后需要加一个基值1
# !现在考虑要将坡地汇流模块单独剖开，只能把产流模块与坡地汇流模块分开
# function ode_sys!(du, u, p, t)
#     itp_input = itp(ceil(t)) # t 2.1->3
#     tmp_u = (itp_input)*(t-u[end]) .* p + u[1:end-1]   # x * dt * p + u
#     tmp_u = circshift(tmp_u, floor(t) - floor(u[end]))
#     tmp_u[end-2] = 0
#     println("u: $u, dt: $(t-u[end]), t-1: $(u[end]), t: $(t)")
#     # println(t-u[end])
#     # println(tmp_u)
#     # println("t: $t, curr u: $u")
#     du[1:end] = vcat(tmp_u, [t-u[end]])
#     # du[end] = t
#     # print("$t" * " : ")
# end

delay = ceil(x4 / delta_t) |> Int
init_states = vcat(zeros(delay)) # states + x0 + t
p = [0.0436344884754979, 0.203199453081548, 0.433360417459522, 0.319805640983432]

prob = DiscreteProblem(ode_sys!, init_states, (0.0, length(input)), p)
sol = solve(prob) # , saveat=1.0
sol.u
f = (x, t) -> (p[1] * sol.u[t][1])


prob = ODEProblem(ode_sys!, init_states, (1.0, length(input)), p)
ode_sol = solve(prob, saveat=1.0) # , saveat=1.0
ode_sol.u

# prob = ODEProblem(ode_sys!, init_states, (1, length(input)), p)
# sol = solve(prob)
# sol.u
# for (idx,v) in enumerate(input)
#     println(p[1]*v+sol.u[idx][1])
# end