"""
Exp-Hydro model
"""
@kwdef mutable struct ExpHydro{T} <: Unit where {T<:Number}
    id::String

    # parameters
    parameters::Dict{Symbol,T}

    # model structure
    elements::Vector{Component}

    # inner variables
    fluxes::ComponentVector{T} = ComponentVector()

    # attribute
    param_names::Vector{Symbol} = [:Tmin, :Tmax, :Df, :Smax, :Qmax, :f]
    input_names::Vector{Symbol} = [:Prcp, :Temp, :Lday]
    output_names::Vector{Symbol} = [:Q]
end

function ExpHydro(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}
    elements = [
        SnowReservoir(id="sr", parameters=parameters, init_states=init_states),
        SoilWaterReservoir(id="wr", parameters=parameters, init_states=init_states)
    ]
    ExpHydro{T}(id=id, elements=elements)
end
