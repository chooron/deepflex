using ModelingToolkit
using DataInterpolations
using DataFrames
using ComponentArrays
using DifferentialEquations
using DataFrames
using CSV

@variables t
const D = Differential(t)

function step_func(x::T; p1=5.0, p2=1.0, p3=0.5) where {T<:Number}
    (tanh(p1 * x) + p2) * p3
end

function data_itp(t, time::AbstractVector, value::AbstractVector)
    itp = LinearInterpolation(value, time, extrapolate=true)
    itp(t)
end

@register_symbolic step_func(x::Num)

function snowfall_func(
    i::(@NamedTuple{prcp::T, temp::T}),
    p::(@NamedTuple{tmin::T}),
    sf::Function) where {T<:Number}
    sf(p.tmin - i.temp) * i.prcp
end

function melt_func(
    i::(@NamedTuple{snw::T, temp::T}),
    p::(@NamedTuple{tmax::T, df::T}),
    sf::Function
)::Union{T,Vector{T}} where {T<:Number}
    snow_water, temp = i[:snw], i[:temp]
    Tmax, Df = p[:tmax], p[:df]
    @.(sf(temp - Tmax) * sf(snow_water) * min(snow_water, Df * (temp - Tmax)))
end

function funcs(; name::Symbol)
    @variables snowfall(t) melt(t) prcp(t) temp(t) temp(t) snw(t) flow(t)
    @parameters tmin tmax df
    eqs = [
        snowfall ~ snowfall_func((prcp=prcp, temp=temp), (tmin=tmin,), step_func)
        melt ~ melt_func((snw=snw, temp=temp), (tmax=tmax, df=df), step_func)
        flow ~ snowfall + melt
    ]
    ODESystem(eqs, t; name=name)
end

function dfuncs(; name::Symbol)
    @variables snw(t) snowfall(t) melt(t)
    @parameters tmin tmax df
    ## SnowWater Reservoir
    eqs = [
        D(snw) ~ snowfall - melt
    ]
    ODESystem(eqs, t; name=name)
end

## start solve
# load data
prcp_data = [3.1, 4.24, 8.02, 15.27, 8.48, 0.0, 0.0, 2.1, 3.55, 0.0]
temp_data = [6.08, 10.53, 11.83, 7.38, 4.8, 5.41, 6.405, 7.675, 5.48, 3.67]

prcp_itp(t) = data_itp(t, collect(1.0:1.0:length(prcp_data)), prcp_data)
temp_itp(t) = data_itp(t, collect(1.0:1.0:length(temp_data)), temp_data)
@register_symbolic prcp_itp(t::Num)
@register_symbolic temp_itp(t::Num)

@variables prcp(t) temp(t)
met = ODESystem([prcp ~ prcp_itp(t)
        temp ~ temp_itp(t)], t, name=:met)

@named sn_funcs = funcs()
@named sn_dfuncs = dfuncs()

conn = compose(ODESystem([
            sn_funcs.prcp ~ met.prcp,
            sn_funcs.temp ~ met.temp,
            sn_dfuncs.melt ~ sn_funcs.melt,
            sn_dfuncs.snowfall ~ sn_funcs.snowfall,
            sn_dfuncs.snw ~ sn_funcs.snw,
        ], t; name=:connnected), met, sn_funcs, sn_dfuncs)

simplified_sys = structural_simplify(conn)

x0 = [simplified_sys.sn_dfuncs.snw => 0.0]
p = [
    simplified_sys.sn_funcs.tmin => -2.092959084,
    simplified_sys.sn_funcs.tmax => 0.175739196,
    simplified_sys.sn_funcs.df => 2.674548848,
]
prob = ODEProblem(simplified_sys, x0, (0.0, 100.0), p)
sol = solve(prob, Tsit5())
