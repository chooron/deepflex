function InfiltrationFlux(
    input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:infiltration;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=infiltration_func,
        smooth_func=smooth_func
    )
end

function infiltration_func(
    i::namedtuple(:snowwater, :liquidwater, :rainfall, :melt),
    p::namedtuple(:whc);
    kw...
)
    sf = get(kw, :smooth_func, step_func)
    @.[sf(i[:liquidwater] - p[:whc] * i[:snowwater]) *
       (i[:rainfall] + i[:melt] + i[:liquidwater] - p[:whc] * i[:snowwater])]
end


function infiltration_func(
    i::namedtuple(:rainfall, :melt),
    p::NamedTuple;
    kw...
)
    @.[i[:rainfall] + i[:melt]]
end

function infiltration_func(
    i::namedtuple(:rainfall),
    p::NamedTuple;
    kw...
)
    [i[:rainfall]]
end

export InfiltrationFlux, infiltration_func