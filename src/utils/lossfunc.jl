function mae(predict::Vector{T}, target::Vector{T}) where T
    sum(abs(target .- predict)) / length(target)
end

function mse(predict::Vector{T}, target::Vector{T}) where T
    sum((target .- predict) .^ T.(2)) / length(target)
end

function rmse(predict::Vector{T}, target::Vector{T}) where T
    sqrt(mse(predict, target))
end

function nse(predict::Vector{T}, target::Vector{T}) where T
    sum((predict .- target) .^ 2) / sum((target .- mean(target)) .^ 2) - 1
end