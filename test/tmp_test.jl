using ComponentArrays
using NamedTupleTools
struct testD <: Any
       parameters::ComponentVector
end

d = testD(ComponentVector(a=1, b=2))
d.parameters.a = 2