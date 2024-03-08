function test1(a::Vector{T}) where {T<:Number}
       for i in a
              if !(i isa AbstractFloat)
                     i = Float32(i)
              end
       end
end

function convert(v::Int32)
       Float32(v)
end

function convert(v::Float32)
       v
end

function test2(a::Vector{T}) where {T<:Number}
       for i in a
              i = convert(i)
       end
end

@btime test1(ones(Int32, 100))
@btime test2(ones(Int32, 100))