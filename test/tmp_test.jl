using ComponentArrays
using ComponentArrays:Axis

function cv_type(arr::List, dtype)
    return ComponentVector{dtype,Vector{dtype},Tuple{Axis{(a=1, b=2)}}}
end

function fun(input::ComponentVector{T,Vector{T},Tuple{Axis{(a=1, b=2)}}}; ) where {T<:Number}
    println(a + b + c)
end

function fun(input::ComponentVector{T,Vector{T},Tuple{Axis{(a=1, c=2)}}}; ) where {T<:Number}
    println(a + b + c)
end

function fun(a; d, e)
    println(a * d * e)
end

tmp_a = 1
tmp_b = 2
tmp_c = 3
tmp_d = 3
tmp_e = 4
fun((tmp_a)...; Dict(:b => tmp_b, :c => tmp_c, :d => tmp_c, :e => tmp_c)...)

@variable $a