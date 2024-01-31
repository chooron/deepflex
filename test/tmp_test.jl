function f1(; kwargs...)
    a = kwargs[:a]
    b = kwargs[:b]
    c = kwargs[:c]
    d = get(kwargs, :d, 4)
    println(a + b + c + d)
end

r = f1(a=1, b=2, c=3)