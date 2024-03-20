function PercolationFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:percolation;
    param_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=percolation_func,
    )
end

function percolation_func(
    i::gen_namedtuple_type([:soilwater], T),
    p::gen_namedtuple_type([:x1], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.((p[:x1]^(-4)) / 4 * ((4 / 9)^(-4)) * (i[:soilwater]^5))
end