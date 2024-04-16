using Optimization
using OptimizationOptimJL
using ComponentArrays
using ComponentArrays:Axis

x0 = [0.0, 0.0]
x_axis = getaxes(ComponentVector(a=1, b=2))

rosenbrock(x, p) = begin
    x_ca = ComponentVector(x, x_axis)
    (p[1] - x_ca.a)^2 + p[2] * (x_ca.b - x_ca.a^2)^2
end

p = [1.0, 100.0]

using OptimizationBBO
prob = OptimizationProblem(rosenbrock, x0, p, lb=ComponentVector([-1.0, -1.0],x_axis), ub=ComponentVector([1.0, 1.0],x_axis))
# prob = OptimizationProblem(rosenbrock, x0, p)
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())

