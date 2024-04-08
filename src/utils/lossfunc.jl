function mae(predict::Vector{T}, target::Vector{T}) where {T<:Number}
    sum(abs(target .- predict)) / length(target)
end

function mse(predict::Vector{T}, target::Vector{T}) where {T<:Number}
    sum((target .- predict) .^ 2) / length(target)
end

function rmse(predict::Vector{T}, target::Vector{T}) where {T<:Number}
    sqrt(mse(predict, target))
end

function nse(predict::Vector{T}, target::Vector{T}) where {T<:Number}
    sum((predict .- target) .^ 2) / sum((target .- mean(target)) .^ 2) - 1
end