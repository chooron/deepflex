"""
Exp-Hydro model
"""
function ExpHydro(; name::String)
    elements = [
        Surface_ExpHydro(name=:sf),
        Soil_ExpHydro(name=:sl),
        Routing_ExpHydro(name=:rt)
    ]
    build_unit(name=name, elements=elements)
end
