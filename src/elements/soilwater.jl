
@kwdef mutable struct SoilWaterReservoir{T<:Number} <: ODEsElement
    id::String

    # flux element
    flux_eles::Vector{ParameterizedElement}()

    # states
    SoilWater::T
    states::Dict{Symbol,Vector{T}} = Dict{Symbol,Vector{T}}()

    # attribute
    input_names::Vector{Symbol} = [:Rain, :Melt, :Temp, :Lday]
    output_names::Vector{Symbol} = [:Pet, :Evap, :Qb, :Qs]
    state_names::Vector{Symbol} = [:SoilWater]
end

function SoilWaterReservoir(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}
    flux_eles = [
        Pet(input_names=[:Temp, :Lday], parameters=Dict()),
        Evap(input_names=[:SoilWater, :Pet], parameters=Dict(:Smax => parameters[:Smax])),
        Baseflow(input_names=[:SoilWater], parameters=Dict(:Smax => parameters[:Smax], :Qmax => parameters[:Qmax], :f => parameters[:f])),
        Surfaceflow(input_names=[:SoilWater], parameters=Dict(:Smax => parameters[:Smax]))
    ]
    SoilWaterReservoir{T}(id=id, SoilWater=init_states[:SoilWater]; flux_eles...)
end

function get_du(ele::SoilWaterReservoir{T}; S::ComponentVector{T}, input::ComponentVector{T}) where {T<:Number}
    fluxes = get_fluxes(ele, S=S, input=input)
    du = @.(fluxes[:Rain] + fluxes[:Melt] - fluxes[:Evap] - fluxes[:Qb] - fluxes[:Qs])
    return ComponentVector(SoilWater=du)
end

function get_fluxes(ele::SoilWaterReservoir{T}; S::ComponentVector{T}, input::ComponentVector{T}) where {T<:Number}
    temp_fluxes = ComponentVector(input; S...)
    for flux_ele in ele.flux_eles
        temp_flux = get_output(flux_ele, input)
        temp_fluxes = ComponentVector(temp_fluxes; temp_flux...)
    end
    return temp_fluxes
end