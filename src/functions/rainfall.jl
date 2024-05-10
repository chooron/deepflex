function RainfallFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:rainfall;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)

    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=rainfall_func,
        smooth_func=smooth_func
    )
end

function rainfall_func(
    i::namedtuple(:prcp, :temp),
    p::namedtuple(:Tmin);
    kw...
)
    sf = get(kw, :smooth_func, step_func)
    @.(sf(i[:temp] - p[:Tmin]) * i[:prcp])
end

function rainfall_func(
    i::namedtuple(:prcp, :pet),
    p::NamedTuple;
    kw...
)
    sf = get(kw, :smooth_func, step_func)
    @.(sf(i[:prcp] - i[:pet]) * (i[:prcp] - i[:pet]))
end

function rainfall_func(
    i::namedtuple(:prcp),
    p::NamedTuple;
    kw...
)
    i[:prcp]
end

function rainfall_func(
    i::namedtuple(:prcp, :temp),
    p::namedtuple(:tt, :tti);
    kw...
)
    tmp_t1 = p[:tt] - 0.5 * p[:tti]
    tmp_t2 = p[:tt] + 0.5 * p[:tti]
    sf = get(kw, :smooth_func, step_func)
    @.(sf(i[:temp] - tmp_t2) * i[:prcp] +
       sf(tmp_t2 - i[:temp]) * sf(i[:temp] - tmp_t1) * i[:prcp] * (i[:temp] - tmp_t1) / p[:tti])
end
