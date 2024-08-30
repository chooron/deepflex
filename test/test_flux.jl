using ComponentArrays
using BenchmarkTools

cv = ComponentVector(a=1, b=2, c=3)
@btime cv[[1,2,3]]

@btime [cv[i] for i in [1,2,3]]