function mae(predict::Vector{<:Number}, target::Vector{<:Number})
    sum(abs(target .- predict)) / length(target)
end

function mse(predict::Vector{<:Number}, target::Vector{<:Number})
    sum((target .- predict) .^ (2)) / length(target)
end

function rmse(predict::Vector{<:Number}, target::Vector{<:Number})
    sqrt(mse(predict, target))
end

function nse(predict::Vector{<:Number}, target::Vector{<:Number})
    sum((predict .- target) .^ 2) / sum((target .- mean(target)) .^ 2) - 1
end