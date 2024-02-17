using ModelingToolkit
using Interpolations
using DifferentialEquations

@variables t
const D = Differential(t)

function step_func(x::T; p1=5.0, p2=1.0, p3=0.5) where {T<:Number}
    (tanh(p1 * x) + p2) * p3
end

@connector function data_itp(input_names::AbstractVector{Symbol}; name::Symbol)
    var_list = [@variables $nm(t) for nm in input_names]
    ODESystem(Equation[], t, [var[1] for var in var_list], []; name=name)
end

function snow_water_reservoir(; name::Symbol)

    ## SnowWater Reservoir
    @variables begin
        t
        prcp(t)
        temp(t)
        lday(t)
        snowfall(t)
        melt(t)
        snw(t)
    end
    @parameters begin
        tmin
        tmax
        df
    end
    D = Differential(t)
    eqs = [
        snowfall ~ step_func(tmin - temp) * prcp,
        melt ~ step_func(temp - tmax) * step_func(snw) * min(snw, df * (temp - tmax)),
        D(snw) ~ snowfall - melt,
    ]
    ODESystem(eqs, t; name=name)
end

function soil_water_reservoir(; name::Symbol)
    ## SnowWater Reservoir
    @variables begin
        t
        prcp(t)
        temp(t)
        lday(t)
        melt(t)
        rainfall(t)
        pet(t)
        evap(t)
        baseflow(t)
        surfaceflow(t)
        flow(t)
        sw(t)
    end
    @parameters begin
        tmin
        qmax
        smax
        f
    end
    D = Differential(t)
    eqs = [
        rainfall ~ step_func(temp - tmin) * prcp,
        pet ~ 29.8 * lday * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2),
        evap ~ step_func(sw) * step_func(sw - smax) * pet + step_func(sw) * step_func(smax - sw) * pet * (sw / smax),
        baseflow ~ step_func(sw) * step_func(sw - smax) * qmax + step_func(sw) * step_func(smax - sw) * qmax * exp(-f * (smax - sw)),
        surfaceflow ~ step_func(sw) * step_func(sw - smax) * (sw - smax),
        flow ~ baseflow + surfaceflow,
        D(sw) ~ rainfall + melt - evap - flow
    ]
    ODESystem(eqs, t; name=name)
end

@named module_1 = snow_water_reservoir()
@named module_2 = soil_water_reservoir()

ODEProblem()