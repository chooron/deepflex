using ModelingToolkit
using Interpolations
using DataFrames
using DifferentialEquations

@variables t
const D = Differential(t)

function step_func(x::T; p1=5.0, p2=1.0, p3=0.5) where {T<:Number}
    (tanh(p1 * x) + p2) * p3
end

@register_symbolic step_func(x)

function snowfall_func(input::(@NamedTuple{prcp::T, temp::T}), paramerter::(@NamedTuple{tmin::T})) where {T<:Number}
    step_func(paramerter.tmin - input.temp) * input.prcp
end

vars = @variables begin
    prcp(t) # precipitation
    temp(t) # temperature
    snowfall(t)
    melt(t)
    snw(t) # snow water
end

ps = @parameters begin
    tmin
    tmax
    df
end

function snow_water_reservoir(; name::Symbol)
    ## SnowWater Reservoir
    eqs = [
        snowfall ~ step_func(tmin - temp) * prcp,
        melt ~ step_func(temp - tmax) * step_func(snw) * min(snw, df * (temp - tmax)),
        D(snw) ~ snowfall - melt,
    ]
    ODESystem(eqs, t; name=name)
end

function rainfall_func(input::(@NamedTuple{prcp::T, temp::T}), paramerter::(@NamedTuple{tmin::T})) where {T<:Number}
    step_func(input.temp - paramerter.tmin) * input.prcp
end

# function soil_water_reservoir(; name::Symbol)
#     ## SnowWater Reservoir
#     @variables prcp(t) [unit = "mm"]
#     @variables begin
#         temp(t)
#         lday(t)
#         melt(t)
#         rainfall(t)
#         pet(t)
#         evap(t)
#         baseflow(t)
#         surfaceflow(t)
#         flow(t)
#         sw(t)
#     end
#     @parameters begin
#         tmin
#         qmax
#         smax
#         f
#     end
#     D = Differential(t)
#     eqs = [
#         rainfall ~ rainfall_func((prcp=prcp, temp=temp), (tmin=tmin,)),
#         # rainfall ~ step_func(temp - tmin) * prcp,
#         pet ~ 29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2),
#         evap ~ step_func(sw) * step_func(sw - smax) * pet + step_func(sw) * step_func(smax - sw) * pet * (sw / smax),
#         baseflow ~ step_func(sw) * step_func(sw - smax) * qmax + step_func(sw) * step_func(smax - sw) * qmax * exp(-f * (smax - sw)),
#         surfaceflow ~ step_func(sw) * step_func(sw - smax) * (sw - smax),
#         flow ~ baseflow + surfaceflow,
#         D(sw) ~ rainfall + melt - evap - flow
#     ]
#     ODESystem(eqs, t; name=name)
# end


@named module_1 = snow_water_reservoir()
simple_module_1 = structural_simplify(module_1)
# @named module_2 = soil_water_reservoir()

## start solve
# import data
using DataFrames
using CSV
using Interpolations

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
data_df = DataFrame(data);
xs = 1:10
flow_vec = data_df[xs, "flow(mm)"]

lday = [0.471, 0.470, 0.46, 0.463, 0.463, 0.46, 0.459, 0.455, 0.452, 0.452]

prcp_data = [3.1, 4.24, 8.02, 15.27, 8.48, 0.0, 0.0, 2.1, 3.55, 0.0]
temp_data = [6.08, 10.53, 11.83, 7.38, 4.8, 5.41, 6.405, 7.675, 5.48, 3.67]

prcp_itp = linear_interpolation(xs, prcp_data)
temp_itp = linear_interpolation(xs, temp_data)

u0 = [
    prcp => prcp_itp,
    temp => temp_itp,
    snw => 0.0,
]
p = [
    tmin => -2.092959084,
    tmax => 0.175739196,
    df => 2.674548848,
]
@named module_1 = snow_water_reservoir()
prob = ODEProblem(module_1, u0, (0.0, 100.0), p)
ModelingToolkit.parameters(module_1)