include("../fluxes/snowfall.jl")
include("../fluxes/rainfall.jl")
include("../fluxes/pet.jl")
include("../fluxes/evap.jl")
include("../fluxes/melt.jl")


@with_kw_noshow struct InterceptionFilter <: ParameterizedElement
    id::String
    Tmin::Union{float,Vector{float}}

    num_upstream::Int = 1
    num_downstream::Int = 1
    input_names::Vector{String} = ["P", "T", "Lday"]
    output_names::Vector{String} = ["Snow", "Rain"]
end

function InterceptionFilter(; id::String, parameters::Dict{String,Union{float,Vector{float}}})
    InterceptionFilter(id, get(parameters, "Tmin", 0.0))
end

function get_output(ele::InterceptionFilter; input::Dict{String,Vector{Number}}, solve::Bool=true)::Dict{String,Vector{Number}}
    flux_snow = snowfall(P=input["P"], T=input["T"], Tmin=ele.Tmin)
    flux_rain = rainfall(P=input["P"], T=input["T"], Tmin=ele.Tmin)
    return [flux_snow, flux_rain]
end

@with_kw_noshow mutable struct SnowReservoir <: ODEsElement
    id::String

    # parameters
    Tmax::float
    Df::float

    # states
    init_states::float
    states::Vector{float}

    # solver
    solver::Any

    # attribute
    num_upstream::Int = 1
    num_downstream::Int = 1
    state_names::Vector{String} = ["SnowWater"]
    input_names::Vector{String} = ["Snow", "Temp"]
    output_names::Vector{String} = ["Melt"]
end

function SnowReservoir(; id::String, parameters::Dict{String,Union{float,Vector{float}}}, init_states::Dict{String,float}, solver::Any)
    SnowReservoir(id, get(parameters, "Tmax", 1.0), get(parameters, "Df", 1.0), get(init_states, "Snow", 0.0), solver)
end

function get_du(ele::SnowReservoir, S::float, input::Dict{String,Number})
    return [Snow - melt(S, get(input, "Temp", 0.0), Tmax=ele.Tmax, Df=ele.Df)]
end

function get_fluxes(ele::SnowReservoir, S::Vector{Number}, input::Dict{String,Vector{Number}})
    return Dict("Melt" => melt(S, get(input, "Temp", 0.0), Tmax=ele.Tmax, Df=ele.Df))
end


@with_kw_noshow mutable struct SoilWaterReservoir <: ODEsElement
    id::String

    # parameters
    Tmin::float
    Smax::float
    Qmax::float
    f::float

    # states
    init_states::float
    states::Vector{float}

    # solver
    solver::Any

    # attribute
    num_upstream::Int = 1
    num_downstream::Int = 1
    state_names::Vector{String} = ["SoilWater"]
    input_names::Vector{String} = ["Rain", "Melt", "Temp", "Lday"]
    output_names::Vector{String} = ["Melt"]

end

function get_du(ele::SnowReservoir, Rain, Melt, T, Lday)
    flux_rain = Rain
    flux_melt = Melt
    pet = pet(T, Lday)
    flux_et = evap(S, pet, ele.Tmax)
    flux_qb = baseflow(S, ele.Smax, ele.Qmax, ele.f)
    flux_qs = surfaceflow(S, ele.Smax)

    du = flux_rain + flux_melt - flux_et - flux_qb - flux_qs
    return [du]
end