@with_kw_noshow mutable struct Unit <: Component
    id::String
    elements::Vector{Element}
    topology::AbstractGraph


    # inner variables
    fluxes::Dict{String,Vector{float}} = nothing
end

