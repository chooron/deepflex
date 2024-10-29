using OrdinaryDiffEq
using DataInterpolations
using Plots
using BenchmarkTools

# 输入流量函数 I(t)
I_vec = Float64[2, 3, 4, 3, 5, 6, 7, 8, 9, 10, 5, 3, 2]
I_vec_copy = copy(I_vec)
timeidx = collect(0:length(I_vec)-1) .* 24
# timeidx = 1:length(I_vec)
# 参数
k = K = 1.0    # 蓄泄常数
x = 0.5    # 权重系数

function route1(I_vec, timeidx, ps)
    itp_func = LinearInterpolation(I_vec, timeidx)
    # 定义 ODE 的右侧函数
    function muskingum_ode!(dS, S, p, t)
        I = itp_func(t)  # 输入流量作为时间的函数
        K, x = p
        S_ = S[1]
        # 计算出流量 Q
        Q = (S_ - K * x * I) / (K * (1 - x))
        # 计算 dS/dt
        dS[1] = I - Q
    end

    # p = (K=K, x=x)
    # 初始条件
    S0 = [K * I_vec[1]]  # 初始储水量

    # 构建 ODE 问题
    prob = ODEProblem(muskingum_ode!, S0, (timeidx[1], timeidx[end]), ps)

    # 求解 ODE 问题
    sol = solve(prob, Tsit5(), saveat=timeidx, dt=24)
    sol_vec = Vector(sol)
    I_vec_2 = [itp_func(t) for t in sol.t]
    Q_vec = @. (sol_vec - K * x * I_vec_2) / (K * (1 - x))
    return Q_vec
end

function route2(I_vec, dt)

    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))

    function msk_prob!(du,u, p, t)
        q_in, q_out = u[1], u[2]
        c0, c1, c2 = p
        new_q = (c0 * I_vec[Int(t)]) + (c1 * q_in) + (c2 * q_out)
        du[1] = I_vec[Int(t)] - q_in
        du[2] = new_q - q_out
    end

    prob = DiscreteProblem(msk_prob!, [I_vec[1], I_vec[1]], (1, length(I_vec)), (c0, c1, c2))
    sol = solve(prob, FunctionMap{true}())
    sol_arr = Array(sol)
    return sol_arr[2, :]
end

function route3(I_vec, dt)
    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))
    q_out = [I_vec[1]]
    for i in 2:length(I_vec)
        q_out = push!(q_out, (c0 * I_vec[i]) + (c1 * I_vec[i-1]) + (c2 * q_out[i-1]))
    end
    return q_out
end

@btime q_out1 = route1(I_vec, timeidx, (K, x))
# q_out2 = route1(I_vec, timeidx, (K, x))
# q_out1 = route1(I_vec, 1:length(I_vec), (K, x))
@btime q_out2 = route2(I_vec, 24)
@btime q_out3 = route3(I_vec, 24)

plot(q_out1, color="red")
plot!(q_out2, color="blue")
plot!(q_out3, color="green")
