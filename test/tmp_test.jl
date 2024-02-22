function func(
    input::Union{@NamedTuple{A::T},@NamedTuple{A::Vector{T}}},
    p::(@NamedTuple{B::T})
) where {T<:Number}
    @info input
    @info p
end

function func2(
    input::Union{Vector{T},T}
) where {T<:Number}
    @info input
end


func((A=[1],), (B=1,))
func2([1])