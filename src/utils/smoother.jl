function ifelse_func(x)
    if x > 0.0
        return 1.0
    else
        return 0.0
    end
end

@register_symbolic ifelse_func(x)

function step_func(x)
    (tanh(5.0 * x) + 1.0) * 0.5
end

@register_symbolic step_func(x)