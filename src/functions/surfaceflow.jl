function SurfaceflowFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:surfaceflow;
    param_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=surfaceflow_func
    )
end

function surfaceflow_func(
    i::namedtuple(:soilwater),
    p::namedtuple(:Smax),
    sf::Function
)
    @.(sf(i[:soilwater]) * sf(i[:soilwater] - p[:Smax]) * (i[:soilwater] - p[:Smax]))
end

function surfaceflow_func(
    i::namedtuple(:surfacerunoff, :prcp),
    p::namedtuple(:Aim),
    sf::Function
)
    @.(i[:surfacerunoff] + p[:Aim] * i[:prcp])
end
