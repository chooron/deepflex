function GR4J(; name::Symbol)
    elements = [
        Surface_GR4J(name=:sf),
        Soil_GR4J(name=:sl),
        Lag_GR4J(name=:lag),
        Routing_GR4J(name=:rt),
    ]
    build_unit(name=name, elements=elements)
end