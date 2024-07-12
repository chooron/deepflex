@doc"""
soilwater movement
Water infiltrates the soil by moving through the surface. Percolation is the movement of water through the soil itself.
Finally, as the water percolates into the deeper layers of the soil, it reaches ground water, which is water below the surface.
The upper surface of this underground water is called the "water table". As you can see in the picture above,
ground water can intersect with surface streams, it can appear at the surface as springs, and it flows generally downhill toward the ocean.
"""

#* infilstraction
function expr(eq::HydroEquation{(:snowwater, :liquidwater, :rainfall, :melt),(:infiltration,),(:whc,)}; kw...)
    snowwater, liquidwater, rainfall, melt = eq.inputs
    whc = first(eq.params)

    [sf(liquidwater - whc * snowwater) *
     (rainfall + melt + liquidwater - whc * snowwater)]
end

function expr(eq::HydroEquation{(:rainfall, :melt),(:infiltration,),()}; kw...)
    rainfall, melt = eq.inputs

    [rainfall + melt]
end

function expr(eq::HydroEquation{(:rainfall,),(:infiltration,),()}; kw...)
    rainfall = first(eq.inputs)

    [rainfall]
end
