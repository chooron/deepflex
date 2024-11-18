@variables a b c d
@parameters k1 k2
flux1 = @hydroflux a ~ k1 * b + k2 * c + d
m = Lux.Chain(Lux.Dense(3, 64), Lux.relu, Lux.Dense(64, 1),  name=:m)
minfo = m([a, b, c])
@nnflux [d] ~ minfo