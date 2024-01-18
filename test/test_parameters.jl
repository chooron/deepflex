using Parameters

@with_kw struct Para3{T<:Number}
    a::T
    b::T
    c::T = a + b
    d::T = a - b
end

Para3{Float64}(a=0.0, b=1.0)