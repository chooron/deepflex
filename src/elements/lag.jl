function Lag_GR4J(; name::Symbol)

    LAGElement(name=name,
        lagtime_dict=Dict(:SlowFlow => :x4, :FastFlow => :x4),
        lag_func=Dict(:SlowFlow => uh_1_half, :FastFlow => uh_2_full)
    )
end

function Lag_HBV(; name::Symbol)

    LAGElement(
        name=name,
        lagtime_dict=Dict(:SlowFlow => :maxbas),
        lag_func=Dict(:Flow => uh_1_half)
    )
end