using DifferentialEquations, RecursiveArrayTools, DiffEqParamEstim
using Optimization, ForwardDiff, OptimizationOptimisers, Optimisers
using DataInterpolations
using DataFrames
using CSV
using DataInterpolations
using SciMLSensitivity
using Zygote
using ComponentArrays
# smoothing step function
step_fct(x) = (tanh(5.0 * x) + 1.0) * 0.5

# snow precipitation
Ps(P, T, Tmin) = step_fct(Tmin - T) * P

# rain precipitation
Pr(P, T, Tmin) = step_fct(T - Tmin) * P

# snow melt
M(S0, T, Df, Tmax) = step_fct(T - Tmax) * step_fct(S0) * minimum([S0, Df * (T - Tmax)])

# evapotranspiration
PET(T, Lday) = 29.8 * Lday * 24.0 * 0.611 * exp((17.3 * T) / (T + 237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1) * step_fct(S1 - Smax) * PET(T, Lday) + step_fct(S1) * step_fct(Smax - S1) * PET(T, Lday) * (S1 / Smax)

# base flow
Qb(S1, f, Smax, Qmax) = step_fct(S1) * step_fct(S1 - Smax) * Qmax + step_fct(S1) * step_fct(Smax - S1) * Qmax * exp(-f * (Smax - S1))

# peak flow
Qs(S1, Smax) = step_fct(S1) * step_fct(S1 - Smax) * (S1 - Smax)
df = DataFrame(CSV.File("data/camels/01013500.csv"))
time = 1:1000
lday_vec = df[time, "dayl(day)"]
prcp_vec = df[time, "prcp(mm/day)"]
temp_vec = df[time, "tmean(C)"]
flow_vec = df[time, "flow(mm)"]
itp_Lday = LinearInterpolation(lday_vec, time)
itp_P = LinearInterpolation(prcp_vec, time)
itp_T = LinearInterpolation(temp_vec, time)

function exp_hydro!(du, u, p, t)
    f, Smax, Qmax, Df, Tmax, Tmin = p
    Lday = itp_Lday(t)
    P = itp_P(t)
    T = itp_T(t)

    Q_out = Qb(u[2], f, Smax, Qmax) + Qs(u[2], Smax)

    du[1] = Ps(P, T, Tmin) - M(u[1], T, Df, Tmax)
    du[2] = Pr(P, T, Tmin) + M(u[1], T, Df, Tmax) - ET(u[2], T, Lday, Smax) - Q_out
end


u0 = [0.0; 1303.004248]
tspan = (1.0, 1000.0)
params = [0.05, 1000.0, 20.0, 3.0, 1.0, -1.0]

prob = ODEProblem(exp_hydro!, u0, tspan, params)
sol = solve(prob, Tsit5(), saveat=1.0, sensealg=ForwardDiffSensitivity())
qout = @.(Qb(sol[2, :], params[1], params[2], params[3]) + Qs(sol[2, :], params[2]))
sum(abs.(qout .- flow_vec))

function loss_func(u, p)
    f, Smax, Qmax, Df, Tmax, Tmin = u
    pt = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
    prob = ODEProblem(exp_hydro!, u0, tspan, u)
    sol = solve(prob, Tsit5(), saveat=1.0)
    Q_out = @.(Qb(sol[2, :], pt.f, pt.Smax, pt.Qmax) + Qs(sol[2, :], pt.Smax))
    sum(abs.(Q_out .- flow_vec))
end

function callback_func(p, l)
    @info l
    false
end

# randomized = VectorOfArray([(sol(t[i]) + 0.01randn(2)) for i in 1:length(t)])
# data = convert(Array, randomized)
# cost_function = build_loss_objective(prob, Tsit5(), L2Loss(t, data),
#     Optimization.AutoForwardDiff(),
#     maxiters=100, verbose=false)
cost_function = Optimization.OptimizationFunction(loss_func, Optimization.AutoZygote())
optprob = Optimization.OptimizationProblem(
    cost_function,
    params,
    # lb=[0.0, 100.0, 10.0, 0.01, 0.0, -3.0],
    # ub=[0.1, 1500.0, 50.0, 5.0, 3.0, 0.0]
)
optsol = solve(optprob, Optimisers.Adam(), maxiters=10, callback=callback_func)