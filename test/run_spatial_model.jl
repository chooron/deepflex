SimpleFlux = HydroModels.SimpleFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
UnitHydroFlux = HydroModels.UnitHydroFlux
NeuralFlux = HydroModels.NeuralFlux
StdMeanNormFlux = HydroModels.StdMeanNormFlux
step_func = HydroModels.step_func

@testset "test grid route hydro model (multiple hydrology nodes based on exp-hydro)" begin
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow flow_routed

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00, s_river=0.0)
    params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09, lag=0.2)
    pas = ComponentVector(initstates=initstates, params=params)

    #! define the snow pack reservoir
    snow_funcs = [
        SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

    #! define the soil water reservoir
    soil_funcs = [
        SimpleFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        SimpleFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
        SimpleFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
    ]
    soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
    soil_ele = HydroBucket(name=:exphydro_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

    #! define the routing method
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    discharge_flux = HydroModels.DischargeRouteFlux(flow, flow_routed)
    discharge_route = HydroModels.GridRoute(name=:exphydro_routed, rfunc=discharge_flux, flwdir=flwdir, positions=positions, subareas=10.0)

    #! define the Exp-Hydro model
    model = HydroModel(name=:exphydro, components=[snow_ele, soil_ele, discharge_route])

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater, :s_river])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :flow_routed])
    @test Set(HydroModels.get_var_names(model)) == Set([:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :snowpack, :soilwater, :s_river, :flow_routed])

    input_ntp_vec = repeat([input_ntp], 9)
    node_names = [Symbol(:node_, i) for i in 1:9]
    node_params = NamedTuple{Tuple(node_names)}(repeat([params], 9))
    node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 9))
    node_pas = ComponentVector(params=node_params, initstates=node_initstates)

    result_mat_vec = model(input_ntp_vec, node_pas, ts, convert_to_ntp=false, ptypes=node_names)
    @test size(result_mat_vec) == (length(HydroModels.get_var_names(model)), length(input_ntp_vec), length(ts))
end

@testset "test vector route hydro model (spatial hydrology model based on vector route and exp-hydro)" begin
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow flow_routed

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00, s_river=0.0)
    params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09, lag=0.2)
    pas = ComponentVector(initstates=initstates, params=params)

    #! define the snow pack reservoir
    snow_funcs = [
        SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

    #! define the soil water reservoir
    soil_funcs = [
        SimpleFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        SimpleFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
        SimpleFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
    ]
    soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
    soil_ele = HydroBucket(name=:exphydro_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

    #! define the routing method
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    discharge_flux = HydroModels.DischargeRouteFlux(flow, flow_routed)
    discharge_route = HydroModels.GridRoute(name=:exphydro_routed, rfunc=discharge_flux, flwdir=flwdir, positions=positions, subareas=10.0)

    #! define the Exp-Hydro model
    model = HydroModel(name=:exphydro, components=[snow_ele, soil_ele, discharge_route])

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater, :s_river])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :flow_routed])
    @test Set(HydroModels.get_var_names(model)) == Set([:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :snowpack, :soilwater, :s_river, :flow_routed])

    input_ntp_vec = repeat([input_ntp], 9)
    node_names = [Symbol(:node_, i) for i in 1:9]
    node_params = NamedTuple{Tuple(node_names)}(repeat([params], 9))
    node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 9))
    node_pas = ComponentVector(params=node_params, initstates=node_initstates)

    result_mat_vec = model(input_ntp_vec, node_pas, ts, convert_to_ntp=false, ptypes=node_names)
    @test size(result_mat_vec) == (length(HydroModels.get_var_names(model)), length(input_ntp_vec), length(ts))
end
