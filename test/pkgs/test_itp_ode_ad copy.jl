using DifferentialEquations # RecursiveArrayTools, DiffEqParamEstim ForwardDiff, 
using Optimization, OptimizationOptimisers
using DataInterpolations
using DataFrames
using CSV
using DataInterpolations
using SciMLSensitivity

using BenchmarkTools
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
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
time = 1:10000
lday_vec = df[time, "dayl(day)"]
prcp_vec = df[time, "prcp(mm/day)"]
temp_vec = df[time, "tmean(C)"]
flow_vec = df[time, "flow(mm)"]
itp_L(t) = LinearInterpolation(lday_vec, time)(t)
itp_P(t) = LinearInterpolation(prcp_vec, time)(t)
itp_T(t) = LinearInterpolation(temp_vec, time)(t)

# eval(Meta.parse("itp_$(var_nm)(t) = LinearInterpolation($(var_nm)_vec, time)(t)"))
# eval(Meta.parse("@register_symbolic itp_$(var_nm)(t)"))

# @register_symbolic itp_L(t)
# @register_symbolic itp_P(t)
# @register_symbolic itp_T(t)

@variables snow(t), soil(t)
@parameters f, Smax, Qmax, Df, Tmax, Tmin

eqs = [
    D(snow) ~ Ps(itp_P(t), itp_T(t), Tmin) - M(snow, itp_T(t), Df, Tmax)
    D(soil) ~ Pr(itp_P(t), itp_T(t), Tmin) + M(snow, itp_T(t), Df, Tmax) - ET(soil, itp_T(t), itp_L(t), Smax) - (Qb(soil, f, Smax, Qmax) + Qs(soil, Smax))
]

u0 = [snow => 0.0, soil => 1303.004248]
@mtkbuild model = ODESystem(eqs, t);
tspan = (time[1], time[end])
params = [f => 0.05, Smax => 1000.0, Qmax => 20.0, Df => 3.0, Tmax => 1.0, Tmin => -1.0]
prob = ODEProblem(model, u0, tspan, params)
sol = solve(prob, Tsit5(); saveat=1.0)