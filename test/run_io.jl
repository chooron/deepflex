ts = collect(1:100)
df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)
pas = ComponentVector(initstates=initstates, params=params)

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
@parameters Tmin Tmax Df Smax f Qmax
@variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow
#! define the snow pack reservoir
snow_funcs = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the soil water reservoir
soil_funcs = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soil_ele = HydroBucket(funcs=soil_funcs, dfuncs=soil_dfuncs)

#! define the Exp-Hydro model
model = HydroModel(name=:exphydro, components=[snow_ele, soil_ele])


@testset "test NamedTupleIOAdapter" begin
    @testset "test single node input" begin
        ioadapter = NamedTupleIOAdapter(model)
        output_ntp = ioadapter(input_ntp, pas)
        @test Set(keys(output_ntp)) == Set(vcat(HydroModels.get_state_names(model), HydroModels.get_output_names(model)))
    end

    @testset "test multiple node input" begin
        ioadapter = NamedTupleIOAdapter(model)
        node_names = [Symbol(:node_, i) for i in 1:10]
        node_params = NamedTuple{Tuple(node_names)}(repeat([params], 10))
        node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 10))
        node_pas = ComponentVector(params=node_params, initstates=node_initstates)
        output_ntp = ioadapter(repeat([input_ntp], 10), node_pas, config=(timeidx=ts, ptypes=node_names))
        @test length(output_ntp) == 10
        @test Set(keys(output_ntp[1])) == Set(vcat(HydroModels.get_state_names(model), HydroModels.get_output_names(model)))
    end
end