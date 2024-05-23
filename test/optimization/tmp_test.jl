using OrdinaryDiffEq, ModelingToolkit
using StructArrays
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using ModelingToolkit: t_nounits as t, D_nounits as D

function lotka_volterra(; name=name)
    states = @variables x(t) = 1.0 y(t) = 1.0
    params = @parameters p1 = 1.5 p2 = 1.0 p3 = 3.0 p4 = 1.0

    eqs = [
        D(x) ~ p1 * x - p2 * x * y,
        D(y) ~ -p3 * y + p4 * x * y
    ]

    return ODESystem(eqs, t, states, params; name=name)
end

@named sys = lotka_volterra()


prob = ODEProblem(structural_simplify(sys), [], (0.0, 10.0), [])
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

using Zygote, SciMLSensitivity

function sum_of_solution(u0, p)
    _prob = remake(prob, u0=[sys.x => u0[1], sys.y => u0[2]], p=p)
    sum(solve(_prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=0.1, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())))
end

u0 = [1.0 1.0]
p = [1.5 1.0 1.0 1.0]
du01, dp1 = Zygote.gradient(sum_of_solution, u0, p)

sa = StructArray((a=rand(1000), b=rand(1000)))
sa2 = StructArray((c=rand(1000), b=rand(1000)))

function test_sa_ntp(s::namedtuple(:a, :b))
    @info s
end

struct A
end

function f1!(v2::Vector, c::Number)
    v2 .*= c
end
a = A()
a1 = [1, 2, 3]

f1!(a1, 3)
a1

function modify!(y::Vector, z::Number)
    # 假设我们要把 y 中的每个元素都加上 z
    for i in eachindex(y)
        y[i] += z
    end
end


function m(av::StructArray, vec::Vector, nm::Symbol)
    getproperty(av, nm) .= vec
end

@btime m(sa, ones(1000), :a)
@btime getproperty(sa, :a)[:] = ones(1000)


axes = getaxes(ComponentVector(a=1, b=2, c=3))

y = StructArray((a=[1, 2, 3, 4], b=[2, 3, 4, 5]))
function f2(x)
    tmp_x = ComponentVector(x, axes)
    y.a .+= tmp_x.a * 2 + tmp_x.b * tmp_x.c
    sum(y.a)
end

Zygote.gradient(f2, [1, 2, 3])