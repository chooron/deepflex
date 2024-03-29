using ComponentArrays
using ModelingToolkit
using DifferentialEquations

include("../../src/DeepFlex.jl")

tmp_input = Float64[1, 2, 3, 4, 5, 3, 2, 1, 4]

rf = DeepFlex.LagFlux(:q, :ql, lag_func=DeepFlex.uh_1_half, param_names=:x4)

for i in tmp_input
    println(rf((q=i,)))
end
# @variables t q(t) ql(t)

# DeepFlex.init_flux!(rf, params=(x4=3.5,))

# # todo ModelingToolkit的注册问题
# # lag_flux(x) = rf((q=x,))
# # @register_symbolic lag_flux(x)
# eqs = [ql ~ DeepFlex.get_output!(rf, (q=q,))]
# println(eqs)
# x1 = DeepFlex.get_output!(rf, (q=q,))
# sys = ODESystem(eqs, t; name=:route)
# equations(sys)
# ! eqs会自动更新，但是sys的不会