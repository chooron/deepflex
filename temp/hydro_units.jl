"""
Exp-Hydro model
"""
function ExpHydro(; name::Symbol)
    elements = [
        Surface_ExpHydro(name=name),
        Soil_ExpHydro(name=name),
        Route_ExpHydro(name=name)
    ]
    HydroUnit(name, elements=elements)
end

function GR4J(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_GR4J(name=name),
        Route_GR4J(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function HBV(; name::Symbol)
    elements = [
        Surface_HBV(name=name),
        Soil_HBV(name=name),
        Route_HBV(name=name)
    ]
    HydroUnit(name, elements=elements)
end

function HyMOD(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_HyMOD(name=name),
        Route_HyMOD(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function XAJ(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_XAJ(name=name),
        Route_XAJ(name=name)
    ]
    HydroUnit(name, elements=elements)
end