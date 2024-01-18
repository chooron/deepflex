@with_kw_noshow struct InterceptionFilter{T<:Number} <: ParameterizedElement
    id::String
    Tmin::T

    num_upstream::Int = 1
    num_downstream::Int = 1
    input_names::Vector{Symbol} = [:Prcp, :Temp]
    output_names::Vector{Symbol} = [:Snow, :Rain]
end

function InterceptionFilter(; id::String, parameters::Dict{Symbol,T}) where {T<:Number}
    InterceptionFilter{T}(id=id, Tmin=get(parameters, :Tmin, 0.0))
end

function get_output(ele::InterceptionFilter; input::Dict{Symbol,Vector{T}}, solve::Bool=true)::Dict{Symbol,Vector{T}} where {T<:Number}
    flux_snow = snowfall(input[:Prcp], input[:Temp], ele.Tmin)
    flux_rain = rainfall(input[:Prcp], input[:Temp], ele.Tmin)
    return Dict(:Snow => flux_snow, :Rain => flux_rain)
end

@with_kw_noshow mutable struct SnowReservoir{T<:Number} <: ODEsElement
    id::String

    # parameters
    Tmax::T
    Df::T

    # states
    init_states::Vector{T}
    states::Vector{T} = []

    # solver
    solver::Any

    # attribute
    num_upstream::Int = 1
    num_downstream::Int = 1
    state_names::Vector{Symbol} = [:SnowWater]
    input_names::Vector{Symbol} = [:Snow, :Temp]
    output_names::Vector{Symbol} = [:Melt]
end

function SnowReservoir(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}, solver::Any) where {T<:Number}
    SnowReservoir{T}(id=id, Tmax=get(parameters, :Tmax, 1.0), Df=get(parameters, :Df, 1.0), init_states=[get(init_states, :SnowWater, 0.0)], solver=solver)
end

function get_du(ele::SnowReservoir, S::Vector{T}, input::Dict{Symbol,T}) where {T<:Number}
    return get(input, :Snow, 0.0) .- melt(S, get(input, :Temp, 0.0), ele.Tmax, ele.Df)
end

function get_fluxes(ele::SnowReservoir, S::Vector{Number}, input::Dict{Symbol,Vector{Number}})
    return Dict(:Melt => melt(S, get(input, :Temp, 0.0), ele.Tmax, ele.Df))
end


@with_kw_noshow mutable struct SoilWaterReservoir{T<:Number} <: ODEsElement
    id::String

    # parameters
    Tmin::T
    Smax::T
    Qmax::T
    f::T

    # states
    init_states::Vector{T}
    states::Vector{T} = []

    # solver
    solver::Any

    # attribute
    num_upstream::Int = 1
    num_downstream::Int = 1
    state_names::Vector{Symbol} = [:SoilWater]
    input_names::Vector{Symbol} = [:Rain, :Melt, :Temp, :Lday]
    output_names::Vector{Symbol} = [:ET, :Qb, :Qs]

end

function SoilWaterReservoir(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}, solver::Any) where {T<:Number}
    println([get(init_states, :SoilWater, 10.0)])
    SoilWaterReservoir{T}(id=id,
        Tmin=get(parameters, :Tmin, 1.0), Smax=get(parameters, :Smax, 1.0),
        Qmax=get(parameters, :Qmax, 0.0), f=get(parameters, :f, 0.0),
        init_states=[get(init_states, :SoilWater, 10.0)], solver=solver)
end

function get_du(ele::SoilWaterReservoir, S::Vector{T}, input::Dict{:Symbol,T}) where {T<:Number}
    flux_rain = input[:Rain]
    flux_melt = input[:Melt]

    pet = pet(input[:Temp], input[:Lday])
    flux_et = evap(S, pet, ele.Tmax)
    flux_qb = baseflow(S, ele.Smax, ele.Qmax, ele.f)
    flux_qs = surfaceflow(S, ele.Smax)

    du = @.flux_rain + flux_melt - flux_et - flux_qb - flux_qs
    return du
end

function get_fluxes(ele::SoilWaterReservoir, S::Vector{T}, input::Dict{Symbol,Vector{T}}) where {T<:Number}
    pet = pet(input[:Temp], input[:Lday])
    flux_et = evap(S, pet, ele.Tmax)
    flux_qb = baseflow(S, ele.Smax, ele.Qmax, ele.f)
    flux_qs = surfaceflow(S, ele.Smax)
    return Dict(:ET => flux_et, :Qb => flux_qb, :Qs => flux_qs)
end

@with_kw_noshow struct FluxAggregator <: BaseElement
    id::String

    num_upstream::Int = 1
    num_downstream::Int = 1
    input_names::Vector{Symbol} = [:Qb, :Qs]
    output_names::Vector{Symbol} = [:Q]
end

function get_output(ele::FluxAggregator, input::Dict{Symbol,Vector{Number}})
    flux_q = input[:Qb] + input[:Qs]
    return Dict(:Q => flux_q)
end

@with_kw mutable struct ExpHydro{T} <: Unit where {T<:Number}
    id::String

    # model structure
    topology::AbstractGraph

    # inner variables
    fluxes::Dict{Symbol,Vector{T}} = Dict()

    # attribute
    input_names::Vector{Symbol} = [:Prcp, :Temp, :Lday]
end

function ExpHydro(; id::String, parameters::Dict{Symbol,T}, init_states::Dict{Symbol,T}) where {T<:Number}

    dag = SimpleDiGraph(4)
    add_edge!(dag, 1, 2)
    add_edge!(dag, 2, 3)
    add_edge!(dag, 3, 4)

    topology = MetaDiGraph(dag)
    set_props!(topology, 1, Dict(:ele => InterceptionFilter(id="ir", parameters=parameters)))
    set_props!(topology, 2, Dict(:ele => SnowReservoir(id="sr", parameters=parameters, init_states=init_states, solver=nothing)))
    set_props!(topology, 3, Dict(:ele => SoilWaterReservoir(id="wr", parameters=parameters, init_states=init_states, solver=nothing)))
    set_props!(topology, 4, Dict(:ele => FluxAggregator(id="fa")))

    ExpHydro{T}(id=id, topology=topology)
end

function get_output(unit::ExpHydro, input::Dict{Symbol,Vector{T}}) where {T<:Number}
    # initialize unit fluxes
    unit.fluxes = input
    # traversal of the directed graph
    for idx in topological_sort(unit.topology)
        tmp_ele = get_prop(unit.topology, idx, :ele)
        tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
        merge!(unit.fluxes, tmp_fluxes)
    end
    return unit.fluxes
end