include("../../src/LumpedHydro.jl")
using Symbolics
using ModelingToolkit
using ModelingToolkit:t_nounits as t
using Lux

@variables temp(t) lday(t) prcp(t) snowpack(t)
@variables pet(t) snowfall(t) melt(t) rainfall(t) infiltration(t)
@parameters Tmin Tmax Df

#* test state flux
snow_funcs = [
    LumpedHydro.SimpleFlux([temp, lday] => [pet],
        flux_exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    LumpedHydro.SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
        flux_exprs=[LumpedHydro.step_func(Tmin - temp) * prcp, LumpedHydro.step_func(temp - Tmin) * prcp]),
    LumpedHydro.SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
        flux_exprs=[LumpedHydro.step_func(temp - Tmax) * LumpedHydro.step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]

funcs_input_ntp = (temp=temp, lday=lday, prcp=prcp,)
funcs_state_ntp = (snowpack=snowpack,)
funcs_params_ntp = (Tmin=Tmin, Tmax=Tmax, Df=Df)
func1,ode_func1 = LumpedHydro.build_ele_func(
    snow_funcs, [LumpedHydro.StateFlux([snowfall] => [melt], snowpack)],
    funcs_input_ntp, funcs_state_ntp, funcs_params_ntp
)
