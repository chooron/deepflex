function RechargeFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:recharge;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=recharge_func,
        smooth_func=smooth_func
    )
end

function recharge_func(
    i::namedtuple(:routingstore),
    p::namedtuple(:x2, :x3, :ω);
    kw...
)
    @.[p[:x2] / (abs(p[:x3])^p[:ω]) * abs(i[:routingstore])^p[:ω]]
end

function recharge_func(
    i::namedtuple(:soilwater, :infiltration),
    p::namedtuple(:fc, :β);
    kw...
)
    @.[(i[:infiltration]) * (i[:soilwater] / p[:fc])^p[:β]]
end

export RechargeFlux, recharge_func