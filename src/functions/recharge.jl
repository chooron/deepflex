function RechargeFlux(input_names::Union{Vector{Symbol},Dict{Symbol,Symbol}},
    output_names::Symbol=:recharge;
    param_names::Vector{Symbol}=Symbol[])
    
    SimpleFlux(
        input_names,
        output_names,
        param_names=param_names,
        func=recharge_func
    )
end

function recharge_func(
    i::gen_namedtuple_type([:routingstore], T),
    p::gen_namedtuple_type([:x2, :x3, :ω], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.(p[:x2] / (p[:x3]^p[:ω]) * i[:routingstore]^p[:ω])
end

function recharge_func(
    i::gen_namedtuple_type([:soilwater, :infiltration], T),
    p::gen_namedtuple_type([:fc, :β], T),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    @.((i[:infiltration]) * (i[:soilwater] / p[:fc])^p[:β])
end