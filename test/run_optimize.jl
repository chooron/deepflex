using Aqua
using CSV
using DataFrames
using Lux
using Test
using ModelingToolkit
using Symbolics
using LuxCore
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
using BenchmarkTools
using Graphs
using Plots
# using HydroModels
include("../src/HydroModels.jl")


# include need add the module name: HydroModels
HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
HydroOptimizer = HydroModels.HydroOptimizer
GradOptimizer = HydroModels.GradOptimizer


@testset "optimize exphydro" begin
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
    # define variables and parameters
    @variables temp lday pet prcp snowfall rainfall snowpack melt
    @parameters Tmin Tmax Df Smax Qmax f

    @variables soilwater pet evap baseflow surfaceflow flow rainfall

    # define model components
    fluxes_1 = [
        HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
    bucket_1 = HydroBucket(name=:surface, funcs=fluxes_1, dfuncs=dfluxes_1)

    fluxes_2 = [
        HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
        HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
    ]
    dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
    bucket_2 = HydroBucket(name=:soil, funcs=fluxes_2, dfuncs=dfluxes_2)
    model = HydroModel(name=:exphydro, components=[bucket_1, bucket_2])

    # load data
    data = CSV.read("data/exphydro/01013500.csv", DataFrame)

    # predefine the parameters
    f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

    # load data
    file_path = "data/exphydro/01013500.csv"
    data = CSV.File(file_path)
    df = DataFrame(data)
    ts = collect(1:10000)
    lday_vec = df[ts, "dayl(day)"]
    prcp_vec = df[ts, "prcp(mm/day)"]
    temp_vec = df[ts, "tmean(C)"]
    flow_vec = df[ts, "flow(mm)"]

    tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
    const_pas = ComponentVector(initstates=ComponentVector(snowpack=0.0, soilwater=1300.0))

    # parameters optimization
    input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
    input_matrix = Matrix(reduce(hcat, collect(input))')
    output = (flow=flow_vec,)
    config = (solver=HydroModels.ODESolver(), interp=LinearInterpolation)

    @testset "HydroOptimizer" begin
        # build optimizer
        hydro_opt = HydroModels.HydroOptimizer(component=model, maxiters=100)
        lb_list = [0.0, 100.0, 10.0, 0.0, 0.0, -3.0]
        ub_list = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
        config = (solver=HydroModels.ODESolver(), interp=LinearInterpolation)
        opt_params, loss_df = hydro_opt([input], [output], tunable_pas=tunable_pas, const_pas=const_pas, config=[config], lb=lb_list, ub=ub_list, return_loss_df=true)
    end

    @testset "GradOptimizer" begin
        config = (solver=HydroModels.ODESolver(sensealg=InterpolatingAdjoint(autodiff=true)), interp=LinearInterpolation)
        grad_opt = HydroModels.GradOptimizer(component=model, maxiters=100)
        opt_params, loss_df = grad_opt([input], [output], tunable_pas=tunable_pas, const_pas=const_pas, config=[config], return_loss_df=true)
    end

    @testset "BatchOptimizer" begin
        config = (solver=HydroModels.ODESolver(sensealg=InterpolatingAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
        batch_opt = HydroModels.BatchOptimizer(component=model, maxiters=100, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote())
        opt_pas, loss_df = batch_opt(
            fill(input, 5), fill(output, 5),
            tunable_pas=tunable_pas,
            const_pas=const_pas, config=fill(config, 5),
            return_loss_df=true
        )
    end
end


@testset "optimize m50" begin
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

    #! parameters in the Exp-Hydro model
    @parameters Tmin Tmax Df Smax f Qmax
    #! parameters in normalize flux
    @parameters snowpack_std snowpack_mean
    @parameters soilwater_std soilwater_mean
    @parameters prcp_std prcp_mean
    @parameters temp_std temp_mean

    #! hydrological flux in the Exp-Hydro model
    @variables prcp temp lday pet rainfall snowfall
    @variables snowpack soilwater lday pet
    @variables melt log_evap_div_lday log_flow
    @variables norm_snw norm_slw norm_temp norm_prcp

    #! define the snow pack reservoir
    snow_funcs = [
        HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

    #! define the ET NN and Q NN
    ep_nn = Lux.Chain(
        Lux.Dense(3 => 16, tanh),
        Lux.Dense(16 => 16, leakyrelu),
        Lux.Dense(16 => 1, leakyrelu),
        name=:epnn
    )
    ep_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), ep_nn))))
    q_nn = Lux.Chain(
        Lux.Dense(2 => 16, tanh),
        Lux.Dense(16 => 16, leakyrelu),
        Lux.Dense(16 => 1, leakyrelu),
        name=:qnn
    )
    q_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))))

    #! get init parameters for each NN
    ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
    q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

    #! define the soil water reservoir
    soil_funcs = [
        #* normalize
        HydroFlux([snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
            [snowpack_mean, soilwater_mean, prcp_mean, temp_mean, snowpack_std, soilwater_std, prcp_std, temp_std],
            exprs=[(var - mean) / std for (var, mean, std) in zip([snowpack, soilwater, prcp, temp],
                [snowpack_mean, soilwater_mean, prcp_mean, temp_mean],
                [snowpack_std, soilwater_std, prcp_std, temp_std]
            )]),
        ep_nn_flux,
        q_nn_flux,
    ]

    state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
    soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr)]
    soil_ele = HydroBucket(name=:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)
    #! define the Exp-Hydro model
    m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele])

    # load data
    data = CSV.read("data/exphydro/01013500.csv", DataFrame)
    # predefine the parameters
    f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
    # load data
    file_path = "data/m50/01013500.csv"
    data = CSV.File(file_path)
    df = DataFrame(data)
    ts = collect(1:10000)
    # cols: Baseflow,Evap,Flow,Infiltration,Lday,Melt,Pet,Prcp,Rainfall,Snowfall,Surfaceflow,Temp,SoilWater,SnowWater
    lday_vec = df[ts, "Lday"]
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    flow_vec = df[ts, "Flow"]

    log_flow_vec = log.(flow_vec)
    log_evap_div_lday_vec = log.(df[ts, "Evap"] ./ lday_vec)
    norm_prcp_vec = (prcp_vec .- mean(prcp_vec)) ./ std(prcp_vec)
    norm_temp_vec = (temp_vec .- mean(temp_vec)) ./ std(temp_vec)
    norm_snw_vec = (df[ts, "SnowWater"] .- mean(df[ts, "SnowWater"])) ./ std(df[ts, "SnowWater"])
    norm_slw_vec = (df[ts, "SoilWater"] .- mean(df[ts, "SoilWater"])) ./ std(df[ts, "SoilWater"])
    nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)

    @testset "ep_nn" begin
        ep_grad_opt = HydroModels.GradOptimizer(component=ep_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
        ep_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(ep_nn_flux)]))')
        ep_output = (log_evap_div_lday=log_evap_div_lday_vec,)

        ep_opt_params, epnn_loss_df = ep_grad_opt(
            [ep_input_matrix], [ep_output],
            tunable_pas=ComponentVector(nn=(epnn=ep_nn_params,)),
            const_pas=ComponentVector(),
            return_loss_df=true
        )
    end

    @testset "q_nn" begin
        q_grad_opt = HydroModels.GradOptimizer(component=q_nn_flux, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=1000)
        q_input_matrix = Matrix(reduce(hcat, collect(nn_input[HydroModels.get_input_names(q_nn_flux)]))')
        q_output = (log_flow=log_flow_vec,)

        q_opt_params, qnn_loss_df = q_grad_opt(
            [q_input_matrix], [q_output],
            tunable_pas=ComponentVector(nn=(qnn=q_nn_params,)),
            const_pas=ComponentVector(),
            return_loss_df=true
        )
    end

    @testset "m50" begin
        m50_opt = HydroModels.GradOptimizer(component=m50_model, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=100)
        config = (solver=HydroModels.ODESolver(sensealg=BacksolveAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
        norm_pas = ComponentVector(
            snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
            snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
        )
        m50_const_pas = ComponentVector(
            initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
            params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...)
        )
        m50_tunable_pas = ComponentVector(
            nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
        )
        m50_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
        m50_opt_params, m50_loss_df = m50_opt(
            [m50_input], [q_output],
            tunable_pas=m50_tunable_pas,
            const_pas=m50_const_pas,
            config=[config],
            return_loss_df=true
        )
    end
end
