function Lag_GR4J(; name::Symbol,
    parameters::ComponentVector{T}) where {T<:Number}

    LAGElement(name=name,
        lag_time=parameters[:x4],
        lag_func=Dict(:SlowFlow => uh_1_half, :FastFlow => uh_2_full)
    )
end

function Lag_HBV(; name::Symbol,
    parameters::ComponentVector{T}) where {T<:Number}

    LAGElement(
        name=name,
        lag_time=parameters[:maxbas],
        lag_func=Dict(:Flow => uh_1_half)
    )
end