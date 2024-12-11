@testset "test hydrobucket" begin
    params = ComponentVector(Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
    init_states = ComponentVector(snowpack=0.0)
    pas = ComponentVector(params=params, initstates=init_states)

    ts = collect(1:1000)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
    dtype = eltype(input[1])
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

    @testset "test hydro element (basic element, Snowpack in Exp-Hydro)" begin
        @variables temp lday prcp pet snowfall rainfall melt snowpack
        @parameters Tmin Tmax Df

        snow_funcs = [
            HydroModels.HydroFlux([temp, lday] => [pet],
                exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
            HydroModels.HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
                exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
            HydroModels.HydroFlux([snowpack, temp] => [melt], [Tmax, Df],
                exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
        ]
        snow_dfuncs = [HydroModels.StateFlux([snowfall] => [melt], snowpack)]
        snow_ele = HydroModels.HydroBucket(funcs=snow_funcs, dfuncs=snow_dfuncs)
        @testset "test hydro element info" begin
            @test Set(HydroModels.get_input_names(snow_ele)) == Set((:temp, :lday, :prcp))
            @test Set(HydroModels.get_param_names(snow_ele)) == Set((:Tmin, :Tmax, :Df))
            @test Set(HydroModels.get_output_names(snow_ele)) == Set((:pet, :snowfall, :rainfall, :melt))
            @test Set(HydroModels.get_state_names(snow_ele)) == Set((:snowpack,))
        end
        config = (timeidx=ts, solver=ManualSolver{true}())
        result = snow_ele(input, pas, config=config)
        ele_state_and_output_names = vcat(HydroModels.get_state_names(snow_ele), HydroModels.get_output_names(snow_ele))
        result = NamedTuple{Tuple(ele_state_and_output_names)}(eachslice(result, dims=1))

        @testset "test first output for hydro element" begin
            snowpack0 = init_states[:snowpack]
            pet0 = snow_funcs[1]([input_ntp.temp[1], input_ntp.lday[1]], ComponentVector(params=ComponentVector()))[1]
            snowfall0, rainfall0 = snow_funcs[2]([input_ntp.prcp[1], input_ntp.temp[1]], ComponentVector(params=(Tmin=params.Tmin,)))
            melt0 = snow_funcs[3]([snowpack0, input_ntp.temp[1]], ComponentVector(params=(Tmax=params.Tmax, Df=params.Df)))[1]
            @test snowpack0 == result.snowpack[1]
            @test snowfall0 == result.snowfall[1]
            @test rainfall0 == result.rainfall[1]
            @test melt0 == result.melt[1]
        end

        @testset "test all of the output" begin
            param_func, nn_param_func = HydroModels._get_parameter_extractors(snow_ele, pas)
            itpfunc_list = map((var) -> LinearInterpolation(var, ts, extrapolate=true), eachrow(input))
            ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]
            du_func = HydroModels._get_du_func(snow_ele, ode_input_func, param_func, nn_param_func)
            solver = ManualSolver{true}()
            initstates_mat = collect(pas[:initstates][HydroModels.get_state_names(snow_ele)])
            #* solve the problem by call the solver
            snowpack_vec = solver(du_func, pas, initstates_mat, ts)[1, :]
            pet_vec = snow_funcs[1](Matrix(reduce(hcat, [input_ntp.temp, input_ntp.lday])'), ComponentVector(params=ComponentVector()))[1, :]
            snow_funcs_2_output = snow_funcs[2](Matrix(reduce(hcat, [input_ntp.prcp, input_ntp.temp])'), ComponentVector(params=(Tmin=params.Tmin,)))
            snowfall_vec, rainfall_vec = snow_funcs_2_output[1, :], snow_funcs_2_output[2, :]
            melt_vec = snow_funcs[3](Matrix(reduce(hcat, [snowpack_vec, input_ntp.temp])'), ComponentVector(params=(Tmax=params.Tmax, Df=params.Df)))[1, :]
            @test reduce(vcat, pet_vec) == collect(result.pet)
            @test reduce(vcat, snowfall_vec) == collect(result.snowfall)
            @test reduce(vcat, rainfall_vec) == collect(result.rainfall)
            @test reduce(vcat, melt_vec) == collect(result.melt)
        end

        @testset "test run with multiple nodes input" begin
            input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input], 10)), (1, 3, 2))
            ndtypes = [Symbol("node_$i") for i in 1:10]
            node_pas = ComponentVector(
                params=NamedTuple{Tuple(ndtypes)}([params for _ in 1:10]),
                initstates=NamedTuple{Tuple(ndtypes)}([init_states for _ in 1:10])
            )
            config = (timeidx=ts, ptypes=ndtypes)
            node_output = snow_ele(input_arr, node_pas, config=config)
            single_output = snow_ele(input, pas, config=config)
            target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
            @test node_output == target_output
        end

        @testset "test build state function" begin
            @variables routingstore exch slowflow_routed routedflow
            @parameters x2, x3
            funcs = [
                HydroModels.HydroFlux([routingstore] => [exch], [x2, x3],
                    exprs=[x2 * abs(routingstore / x3)^3.5]),
                HydroModels.HydroFlux([routingstore, slowflow_routed, exch] => [routedflow], [x3],
                    exprs=[x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5]),
            ]
            dfunc = HydroModels.StateFlux([slowflow_routed, exch] => [routedflow], routingstore)
            hydro_bucket = HydroModels.HydroBucket(funcs=funcs, dfuncs=[dfunc])
            # meta = HydroModels.HydroMeta(name=:test, inputs=[:slowflow_routed], states=[:routingstore])
            flux_func, state_func = HydroModels.build_ele_func(funcs, [dfunc], hydro_bucket.meta)
            rgt, slg = 10.0, 20.0
            exch = 2.42 * abs(rgt / 69.63)^3.5
            routedflow = 69.63^(-4) / 4 * (rgt + slg + exch)^5
            @test flux_func([slg, rgt], [2.42, 69.63], [], 1) ≈ [exch, routedflow]
            @test state_func([slg], [rgt], [2.42, 69.63], [], 1) ≈ [exch + slg - 69.63^(-4) / 4 * (rgt + slg + exch)^5]
        end
    end
end