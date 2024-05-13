function PercolationFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:percolation;
    param_names::Vector{Symbol}=Symbol[],
    smooth_func::Function=step_func)
    
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=percolation_func,
        smooth_func=smooth_func
    )
end

function percolation_func(
    i::namedtuple(:soilwater),
    p::namedtuple(:x1);
    kw...
)
    @.((p[:x1]^(-4)) / 4 * ((4 / 9)^(4)) * (i[:soilwater]^5))
end

export PercolationFlux, percolation_func