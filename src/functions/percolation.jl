#* percolation
function expr(eq::HydroEquation{(:soilwater,),(:percolation,),(:x1,)}; kw...)
    soilwater = first(eq.inputs)
    x1 = first(eq.params)

    [(abs(x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]
end