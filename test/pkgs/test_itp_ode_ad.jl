using OrdinaryDiffEq
using Optimization
using DataInterpolations
using DataFrames
using CSV
using DataInterpolations
using SciMLSensitivity
using Zygote
using ComponentArrays
using BenchmarkTools
using Symbolics
using ModelingToolkit: t_nounits as t
using ModelingToolkit: D_nounits as D
using ModelingToolkit
using OptimizationOptimisers
using NamedTupleTools

function setup_prob(
    sys::ODESystem,
    prob::ODEProblem;
    pas::ComponentVector,
)
    params, init_states = pas[:params],pas[:initstates]
    #* setup init states
    u0 = [getproperty(sys, nm) => init_states[nm] for nm in keys(init_states)]
    u0var = ModelingToolkit.varmap_to_vars(u0, unknowns(sys))
    #* setup parameters
    p = [getproperty(sys, Symbol(nm)) => params[Symbol(nm)] for nm in ModelingToolkit.parameters(base_sys)]
    # pvar = ModelingToolkit.MTKParameters(sys, p, u0var)
    new_prob = remake(prob, p=p, u0=u0var)
    new_prob
end

function merge_ca(ca1::ComponentArray, ca2::ComponentArray, key::Symbol=:param)
    #* used for merge data
    share_keys = intersect(keys(ca1), keys(ca2))
    new_ca1 = ComponentVector{promote_type(eltype(ca1),eltype(ca2))}() # {eltype(ca1)}
    if length(share_keys) > 0
        for key in share_keys
            if !(typeof(ca1[key]) <: ComponentArray) | !(typeof(ca2[key]) <: ComponentArray)
                new_ca1 = ComponentVector(new_ca1; ComponentVector([key=>ca2[key]])...)
            else
                merged_ca = merge_ca(ca1[key], ca2[key], key)
                new_ca1 = ComponentVector(new_ca1; merged_ca...)
            end
        end
    end
    for key in filter(k -> !(k in share_keys), keys(ca1))
        new_ca1 = ComponentVector(new_ca1; ca1[[key]]...)
    end
    for key in filter(k -> !(k in share_keys), keys(ca2))
        new_ca1 = ComponentVector(new_ca1; ca2[[key]]...)
    end
    ComponentVector(namedtuple([key],[new_ca1]))
end

# smoothing step function
step_fct(x) = ((tanh(5.0 * x) + 1.0) * 0.5)

# snow precipitation
Ps(P, T, Tmin) = (step_fct(Tmin - T) * P)

# rain precipitation
Pr(P, T, Tmin) = (step_fct(T - Tmin) * P)

# snow melt
M(S0, T, Df, Tmax) = (step_fct(T - Tmax) * step_fct(S0) * minimum([S0, Df * (T - Tmax)]))

# evapotranspiration
PET(T, Lday) = (29.8 * Lday * 24.0 * 0.611 * exp((17.3 * T) / (T + 237.3)) / (T + 273.2))
ET(S1, T, Lday, Smax) = (step_fct(S1) * step_fct(S1 - Smax) * PET(T, Lday) + step_fct(S1) * step_fct(Smax - S1) * PET(T, Lday) * (S1 / Smax))

# base flow
Qb(S1, f, Smax, Qmax) = (step_fct(S1) * step_fct(S1 - Smax) * Qmax + step_fct(S1) * step_fct(Smax - S1) * Qmax * exp(-f * (Smax - S1)))

# peak flow
Qs(S1, Smax) = (step_fct(S1) * step_fct(S1 - Smax) * (S1 - Smax))

time = 1:1000
df = DataFrame(CSV.File("data/camels/01013500.csv"))
lday_vec = df[time, "dayl(day)"]
prcp_vec = df[time, "prcp(mm/day)"]
temp_vec = df[time, "tmean(C)"]
flow_vec = df[time, "flow(mm)"]

itp_Lday(t) = LinearInterpolation(lday_vec, time)(t)
itp_P(t) = LinearInterpolation(prcp_vec, time)(t)
itp_T(t) = LinearInterpolation(temp_vec, time)(t)

@register_symbolic itp_Lday(t)
@register_symbolic itp_P(t)
@register_symbolic itp_T(t)

@parameters f, Smax, Qmax, Df, Tmax, Tmin
@variables snw(t), slw(t)
eqs = [
    D(snw) ~ Ps(itp_P(t), itp_T(t), Tmin) - M(snw, itp_T(t), Df, Tmax),
    D(slw) ~ Pr(itp_P(t), itp_T(t), Tmin) + M(snw, itp_T(t), Df, Tmax) - ET(slw, itp_T(t), itp_Lday(t), Smax) - (Qb(slw, f, Smax, Qmax) + Qs(slw, Smax))
]
base_sys = ODESystem(eqs, t, name=:sys)
sys = structural_simplify(base_sys)
params_list = [0.05, 100.0, 20.0, 3.0, 1.0, -1.0]
u0_list = [snw => 0.0, slw => 1000.0]
prob = ODEProblem(sys, u0_list, (1.0, 100.0), params_list)
input_ca = ComponentVector(lday=lday_vec, prcp=prcp_vec)
p_axes = getaxes(ComponentVector(params=(f=0.05, Smax=100.0, Qmax=20.0, Df=3.0, Tmax=1.0, Tmin=-1.0)))
init_states = ComponentVector(initstates=(snw=0.0, slw=1000.0))
sol = solve(prob, BS3(), dt=1.0)

function loss_func(u, p)
    params = ComponentVector(u, p_axes)
    pas = merge_ca(params, init_states)[:param]
    new_prob = setup_prob(sys, prob, pas=pas)
    sol = solve(new_prob, BS3(), saveat=1.0)
    Q_out = @.(Qb(sol[2, :], u[1], u[2], u[3]) + Qs(sol[2, :], u[3]))
    result = ComponentVector(input_ca; ComponentVector(q=Q_out)...)
    sum(result[:q])
end

function callback_func(p, l)
    @info l
    false
end


cost_function = Optimization.OptimizationFunction(loss_func, Optimization.AutoZygote()) # AutoReverseDiff, AutoZygote
optprob = Optimization.OptimizationProblem(
    cost_function,
    params_list,
)
sol = solve(optprob, Adam(1e-2), maxiters=10, callback=callback_func)