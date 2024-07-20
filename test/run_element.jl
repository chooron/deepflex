params = ComponentVector(Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
init_states = ComponentVector(snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states)

ts = collect(1:100)
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"]);
dtype = eltype(input[1]);

@testset "test hydro element (basic element, Snowpack in Exp-Hydro)" begin
    @variables temp(t) lday(t) prcp(t) pet(t) snowfall(t) rainfall(t) melt(t) snowpack(t)
    @parameters Tmin Tmax Df

    snow_funcs = [
        SimpleFlux([temp, lday] => [pet],
            flux_exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
            flux_exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
            flux_exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroElement(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)
    @testset "test hydro element info" begin
        @test Set(LumpedHydro.get_input_names(snow_ele)) == Set((:temp, :lday, :prcp))
        @test Set(LumpedHydro.get_param_names(snow_ele)) == Set((:Tmin, :Tmax, :Df))
        @test Set(LumpedHydro.get_output_names(snow_ele)) == Set((:pet, :snowfall, :rainfall, :melt))
        @test Set(LumpedHydro.get_state_names(snow_ele)) == Set((:snowpack,))
    end

    result = snow_ele(input, pas)

    @testset "test first output for hydro element" begin
        snowpack0 = init_states[:snowpack]
        pet0 = snow_funcs[1]([input.temp[1] input.lday[1]], dtype[])[1]
        snowfall0, rainfall0 = snow_funcs[2]([input.prcp[1] input.temp[1]], [params.Tmin])
        melt0 = snow_funcs[3]([snowpack0 input.temp[1]], [params.Tmax, params.Df])[1]
        @test snowpack0 == result.snowpack[1]
        @test snowfall0 == result.snowfall[1]
        @test rainfall0 == result.rainfall[1]
        @test melt0 == result.melt[1]
    end

    # todo
    @testset "test build state function" begin
        
    end

    # todo
    @testset "test modify element" begin
        
    end

    @testset "test ode solved results" begin
        prcp_itp = LinearInterpolation(input.prcp, ts)
        temp_itp = LinearInterpolation(input.temp, ts)

        function snowpack_bucket!(du, u, p, t)
            snowpack_ = u[1]
            Df, Tmax, Tmin = p
            prcp_, temp_ = prcp_itp(t), temp_itp(t)
            snowfall_ = step_func(Tmin - temp_) * prcp_
            melt_ = step_func(temp_ - Tmax) * step_func(snowpack_) * min(snowpack_, Df * (temp_ - Tmax))
            du[1] = snowfall_ - melt_
        end
        prob = ODEProblem(snowpack_bucket!, [init_states.snowpack], (ts[1], ts[end]), collect(params))
        sol = solve(prob, Tsit5(), saveat=ts, reltol=1e-3, abstol=1e-3)
        num_u = length(prob.u0)
        manual_result = [sol[i, :] for i in 1:num_u]
        pkg_result = LumpedHydro.solve_prob(snow_ele, input=input, params=params, init_states=init_states)
        @test manual_result == pkg_result
    end

    @testset "test all of the output" begin
        snowpack_vec = LumpedHydro.solve_prob(snow_ele, input=input, params=params, init_states=init_states)[1]
        pet_vec = snow_funcs[1](reduce(hcat, [input.temp, input.lday]), dtype[])[:, 1]
        snow_funcs_2_output = snow_funcs[2](reduce(hcat, [input.prcp, input.temp]), [params.Tmin])
        snowfall_vec, rainfall_vec = snow_funcs_2_output[:, 1], snow_funcs_2_output[:, 2]
        melt_vec = snow_funcs[3](reduce(hcat, [snowpack_vec, input.temp]), [params.Tmax, params.Df])[:, 1]
        @test pet_vec == collect(result.pet)
        @test snowfall_vec == collect(result.snowfall)
        @test rainfall_vec == collect(result.rainfall)
        @test melt_vec == collect(result.melt)
    end
end