function RenameFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
)
    flux_exprs = [var for var in fluxes[1]]
    SimpleFlux(fluxes, Num[], exprs=flux_exprs)
end