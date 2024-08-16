params = ComponentVector(Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
init_states = ComponentVector(snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states)

ts = collect(1:100)
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"]);
input = reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))'
dtype = eltype(input[1]);

@testset "test hydro element (basic element, Snowpack in Exp-Hydro)" begin
    @variables temp lday prcp pet snowfall rainfall melt snowpack
    @parameters Tmin Tmax Df

    snow_funcs = [
        SimpleFlux([temp, lday] => [pet],
            exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
            exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df],
            exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)
    @testset "test hydro element info" begin
        @test Set(LumpedHydro.get_input_names(snow_ele)) == Set((:temp, :lday, :prcp))
        @test Set(LumpedHydro.get_param_names(snow_ele)) == Set((:Tmin, :Tmax, :Df))
        @test Set(LumpedHydro.get_output_names(snow_ele)) == Set((:pet, :snowfall, :rainfall, :melt))
        @test Set(LumpedHydro.get_state_names(snow_ele)) == Set((:snowpack,))
    end

    result = snow_ele(input, pas, timeidx=ts)
    ele_state_and_output_names = vcat(LumpedHydro.get_state_names(snow_ele), LumpedHydro.get_output_names(snow_ele))
    result = NamedTuple{Tuple(ele_state_and_output_names)}(eachslice(result, dims=1))
    @testset "test first output for hydro element" begin
        snowpack0 = init_states[:snowpack]
        pet0 = snow_funcs[1]([input_ntp.temp[1], input_ntp.lday[1]], dtype[])[1]
        snowfall0, rainfall0 = snow_funcs[2]([input_ntp.prcp[1], input_ntp.temp[1]], [params.Tmin])
        melt0 = snow_funcs[3]([snowpack0, input_ntp.temp[1]], [params.Tmax, params.Df])[1]
        @test snowpack0 == result.snowpack[1]
        @test snowfall0 == result.snowfall[1]
        @test rainfall0 == result.rainfall[1]
        @test melt0 == result.melt[1]
    end


    # # todo
    # @testset "test modify element" begin

    # end

    @testset "test ode solved results" begin
        prcp_itp = LinearInterpolation(input_ntp.prcp, ts)
        temp_itp = LinearInterpolation(input_ntp.temp, ts)

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
        pkg_result = LumpedHydro.solve_single_prob(snow_ele, input=input, pas=pas, timeidx=ts)
        @test manual_result[1] == pkg_result[1, :]
    end

    @testset "test all of the output" begin
        snowpack_vec = LumpedHydro.solve_single_prob(snow_ele, input=input, pas=pas, timeidx=ts)[1, :]
        pet_vec = reduce(hcat, snow_funcs[1].(eachslice(reduce(hcat, [input_ntp.temp, input_ntp.lday]), dims=1), Ref(dtype[])))[1,:]
        snow_funcs_2_output = reduce(hcat, snow_funcs[2].(eachslice(reduce(hcat, [input_ntp.prcp, input_ntp.temp]), dims=1), Ref([params.Tmin])))
        snowfall_vec, rainfall_vec = snow_funcs_2_output[1, :], snow_funcs_2_output[2, :]
        melt_vec = reduce(hcat, snow_funcs[3].(eachslice(reduce(hcat, [snowpack_vec, input_ntp.temp]), dims=1), Ref([params.Tmax, params.Df])))[1, :]
        @test reduce(vcat, pet_vec) == collect(result.pet)
        @test reduce(vcat, snowfall_vec) == collect(result.snowfall)
        @test reduce(vcat, rainfall_vec) == collect(result.rainfall)
        @test reduce(vcat, melt_vec) == collect(result.melt)
    end

    # #todo
    # @testset "test build bucket function" begin end
end