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
    i::gen_namedtuple_type([:soilwater], T),
    p::gen_namedtuple_type([:Smax], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(sf(i[:soilwater]) * sf(i[:soilwater] - p[:Smax]) * (i[:soilwater] - p[:Smax]))
end

function surfaceflow_func(
    i::gen_namedtuple_type([:surfacerunoff, :prcp], T),
    p::gen_namedtuple_type([:Aim], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(i[:surfacerunoff] + p[:Aim] * i[:prcp])
end
