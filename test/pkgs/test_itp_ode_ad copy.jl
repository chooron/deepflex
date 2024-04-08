using DifferentialEquations # RecursiveArrayTools, DiffEqParamEstim ForwardDiff, 
using Optimization, OptimizationOptimisers
using DataInterpolations
using DataFrames
using CSV
using DataInterpolations
using SciMLSensitivity
# using Zygote
# using ComponentArrays , Optimisers
using BenchmarkTools
using ModelingToolkit
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
itp_Lday(t) = LinearInterpolation(lday_vec, time)(t)
itp_P(t) = LinearInterpolation(prcp_vec, time)(t)
itp_T(t) = LinearInterpolation(temp_vec, time)(t)

@variables t, snow(t), soil(t)
@parameters f, Smax, Qmax, Df, Tmax, Tmin
@register_symbolic itp_Lday(t)
@register_symbolic itp_P(t)
@register_symbolic itp_T(t)

eqs = [
    snow ~ Ps(itp_P(t), itp_T(t), Tmin) - M(snow, itp_T(t), Df, Tmax),
    soil ~ Pr(itp_P(t), itp_T(t), Tmin) + M(snow, itp_T(t), Df, Tmax) - ET(soil, itp_T(t), itp_Lday(t), Smax) - (Qb(soil, f, Smax, Qmax) + Qs(soil, Smax))
]

u0 = [snow => 0.0, soil => 1303.004248]
model = structural_simplify(ODESystem(eqs, t, name=:temp));
tspan = (time[1], time[end])
# params = [var => value for (var, value) in zip(parameters(model), 1:length(parameters(model)))]
params = [f => 0.05, Smax => 1000.0, Qmax => 20.0, Df => 3.0, Tmax => 1.0, Tmin => -1.0]
# params = [var => value for (var, value) in zip(parameters(model), 1:length(parameters(model)))]
idxs = ModelingToolkit.varmap_to_vars(params, parameters(model))
prob = samefile(model, u0, tspan, params)
sol = solve(prob, FBDF(), saveat=1.0)

# function loss_func(u, p)
#     u0 = [snow => 0.0, soil => 1303.004248]
#     model = structural_simplify(ODESystem(eqs, t, name=:temp))
#     tspan = (time[1], time[end])
#     # params = [var => value for (var, value) in zip(parameters(model), 1:length(parameters(model)))]
#     params = [f => 0.05, Smax => 1000.0, Qmax => 20.0, Df => 3.0, Tmax => 1.0, Tmin => -1.0]
#     # params = [var => value for (var, value) in zip(parameters(model), 1:length(parameters(model)))]
#     idxs = ModelingToolkit.varmap_to_vars(params, parameters(model))
#     prob = ODEProblem(model, u0, tspan, params)
#     sol = solve(prob, FBDF(), saveat=1.0, p=u, abstol=1e-3, reltol=1e-3)
#     Q_out = @.(Qb(sol[2, :], u[1], u[2], u[3]) + Qs(sol[2, :], u[2]))
#     loss = sum(abs.(Q_out .- flow_vec))
#     println(loss)
#     loss
# end

# cost_function = Optimization.OptimizationFunction(loss_func, Optimization.AutoModelingToolkit())
# optprob = Optimization.OptimizationProblem(
#     cost_function,
#     [3.0, 1.0, -1.0, 1000.0, 0.05, 20.0],
# )
# sol = solve(optprob, Optimisers.Adam(1e-2), maxiters=10)

# qout = @.(Qb(sol[2, :], params[1], params[2], params[3]) + Qs(sol[2, :], params[2]))
# sum(abs.(qout .- flow_vec))

# function loss_func(u, p)
#     f, Smax, Qmax, Df, Tmax, Tmin = u
#     pt = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
#     prob = ODEProblem(exp_hydro!, u0, tspan, u)
#     sol = solve(prob, BS3(), saveat=1.0)
#     Q_out = @.(Qb(sol[2, :], pt.f, pt.Smax, pt.Qmax) + Qs(sol[2, :], pt.Smax))
#     loss = sum(abs.(Q_out .- flow_vec))
#     loss
# end

# function callback_func(p, l)
#     @info l
#     false
# end

# # randomized = VectorOfArray([(sol(t[i]) + 0.01randn(2)) for i in 1:length(t)])
# # data = convert(Array, randomized)
# # cost_function = build_loss_objective(prob, Tsit5(), L2Loss(t, data),
# #     Optimization.AutoForwardDiff(),
# #     maxiters=100, verbose=false)
# cost_function = Optimization.OptimizationFunction(loss_func, Optimization.AutoForwardDiff())
# optprob = Optimization.OptimizationProblem(
#     cost_function,
#     params,
#     # lb=[0.0, 100.0, 10.0, 0.01, 0.0, -3.0],
#     # ub=[0.1, 1500.0, 50.0, 5.0, 3.0, 0.0]
# )
# @btime sol = solve(optprob, Optimisers.Adam(1e-2), maxiters=10, callback=callback_func)