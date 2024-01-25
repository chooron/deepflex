"""

"""
@kwdef mutable struct ExpHydro{T} <: Unit where {T<:Number}
    id::String

    # model structure
    structure::AbstractGraph

    # inner variables
    fluxes::ComponentVector{T}=ComponentVector()

    # attribute
    param_names::Vector{Symbol} = [:Tmin, :Tmax, :Df, :Smax, :Qmax, :f]
    input_names::Vector{Symbol} = [:Prcp, :Temp, :Lday]
    output_names::Vector{Symbol} = [:Q]
end

function ExpHydro(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}

    dag = SimpleDiGraph(4)
    add_edge!(dag, 1, 2)
    add_edge!(dag, 2, 3)
    add_edge!(dag, 3, 4)

    structure = MetaDiGraph(dag)
    set_props!(structure, 1, Dict(:ele => InterceptionFilter(id="ir", parameters=parameters)))
    set_props!(structure, 2, Dict(:ele => SnowReservoir(id="sr", parameters=parameters, init_states=init_states, solver=nothing)))
    set_props!(structure, 3, Dict(:ele => SoilWaterReservoir(id="wr", parameters=parameters, init_states=init_states, solver=nothing)))
    set_props!(structure, 4, Dict(:ele => FluxAggregator(id="fa")))

    ExpHydro{T}(id=id, structure=structure)
end
