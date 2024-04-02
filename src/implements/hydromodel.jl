"""
Exp-Hydro model
"""
function ExpHydro(; name::Symbol)
    elements = [
        Surface_ExpHydro(name=name),
        Soil_ExpHydro(name=name),
        Slope_ExpHydro(name=name)
    ]
    HydroUnit(name, elements=elements)
end

function GR4J(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_GR4J(name=name),
        Slope_GR4J(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function HBV(; name::Symbol)
    elements = [
        Surface_HBV(name=name),
        Soil_HBV(name=name),
        Slope_HBV(name=name)
    ]
    HydroUnit(name, elements=elements)
end

function HyMOD(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_HyMOD(name=name),
        Slope_HyMOD(name=name),
    ]
    HydroUnit(name, elements=elements)
end

function XAJ(; name::Symbol)
    elements = [
        Surface_GR4J(name=name),
        Soil_XAJ(name=name),
        Slope_XAJ(name=name)
    ]
    HydroUnit(name, elements=elements)
end