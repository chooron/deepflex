"""
Implement for [Improving hydrologic models for predictions and process understanding using neural ODEs](https://hess.copernicus.org/articles/26/5085/2022/)
"""

@kwdef mutable struct M50{T} <: Unit where {T<:Number}
    id::String

    # model structure
    structure::AbstractGraph

    # inner variables
    fluxes::ComponentVector{T}=ComponentVector()

    # attribute
    param_names::Vector{Symbol} = [:Tmin, :Tmax, :Df, :Smax, :Qmax, :f]
    input_names::Vector{Symbol} = [:Prcp, :Temp, :Lday, :SnowWater, :SoilWater]
    output_names::Vector{Symbol} = [:Q]
end

function M50(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}
    elements = [
        SnowReservoir(id="sr", parameters=parameters, init_states=init_states),
        SoilWaterReservoir(id="wr", parameters=parameters, init_states=init_states)
    ]
    ExpHydro{T}(id=id, elements=elements)
end

