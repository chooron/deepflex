# include("../../../src/HydroModels.jl")


using Random, KolmogorovArnold,ComponentArrays
using Lux
rng = Random.default_rng()

in_dim, out_dim, grid_len = 4, 4, 8
layer = KDense(in_dim, out_dim, grid_len)
p, st = Lux.setup(rng, layer)

length(collect(ComponentArray(p)))

x = rand32(rng, in_dim, 10)
y = layer(x, p, st)


