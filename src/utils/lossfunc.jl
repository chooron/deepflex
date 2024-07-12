function mae(predict::AbstractArray, target::AbstractArray)
    sum(abs(target .- predict)) / length(target)
end

function mse(predict::AbstractArray, target::AbstractArray)
    sum((target .- predict) .^ (2)) / length(target)
end

function rmse(predict::AbstractArray, target::AbstractArray)
    sqrt(mse(predict, target))
end

function nse(predict::AbstractArray, target::AbstractArray)
    sum((predict .- target) .^ 2) / sum((target .- mean(target)) .^ 2)
end