using ComponentArrays

include("../src/DeepFlex.jl")

x4 = 3.5
tmp_input = ComponentVector(Q9=Float64[1, 2, 3, 4, 5, 3, 2, 1, 4])

rf = DeepFlex.RoutingFlux([:Q9], lag_time=x4, lag_func=DeepFlex.unit_hydro)

@time tmp_output = rf(tmp_input)
