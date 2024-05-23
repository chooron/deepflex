function MeltFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:melt;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=melt_func,
        smooth_func=smooth_func
    )
end

function melt_func(
    i::namedtuple(:snowwater, :temp),
    p::namedtuple(:Tmax, :Df);
    kw...
)
    sf = get(kw, :smooth_func, step_func)
    @.[sf(i[:temp] - p[:Tmax]) * sf(i[:snowwater]) * min(i[:snowwater], p[:Df] * (i[:temp] - p[:Tmax]))]
end

function melt_func(
    i::namedtuple(:temp),
    p::namedtuple(:cfmax, :ttm);
    kw...
)
    sf = get(kw, :smooth_func, step_func)
    @.[sf(i[:temp] - p[:ttm]) * (i[:temp] - p[:ttm]) * p[:cfmax]]
end

export MeltFlux, melt_func