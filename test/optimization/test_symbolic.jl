using Symbolics
using ModelingToolkit

@variables t a(t) b(t)
@parameters p1 p2

eqs = [a * p1 + b * p2]
bf = build_function(eqs[1], [a, b], [p1, p2], expression=Val{false})
bf.([[1, 2], [3, 2], [3, 4]], Ref([2, 4]))

bf.([(a=1, b=2), (b=3, a=4), (a=5, b=6)], Ref([2, 4]))


@variables t lday(t) temp(t)
eq1 = (29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2))
bf = build_function(eq1, [lday, temp], expression=Val{false})

v1 = [[1, 2, 4, 3], [2, 3, 4, 3], [3, 2, 1, 3]]

function pront(v)
    println(v)
    println(typeof(vec(v).*2.0))
end

pront.(eachslice(hcat(v1...), dims=1))

vs = eachslice(hcat(v1...), dims=1)