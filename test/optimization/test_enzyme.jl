using Zygote
using ComponentArrays
using StructArrays


axes = getaxes(ComponentVector(a=1, b=2, c=3))

function f2(x)
    y = (a=[1, 2, 3, 4], b=[2, 3, 4, 5])
    tmp_x = ComponentVector(x, axes)
    dt = tmp_x.a * 2 .+ tmp_x.b * tmp_x.c
    Base.setindex
    y = merge(y, (a=y.a .+ dt,))
    sum(y.a)
end

Zygote.gradient(f2, [1, 2, 3])