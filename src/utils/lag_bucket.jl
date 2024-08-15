
struct LagBucket <: AbstractLagBucket
    "the name of hydrological computation bucket "
    name::Symbol
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple
    "lag hydrological fluxes, used to flood routing"
    lfuncs::Vector

    function LagBucket(
        name::Symbol;
        lfuncs::Vector{<:AbstractLagFlux},
    )
        #* Extract all variable names of funcs and dfuncs
        ele_input_names, ele_output_names = reduce(union, get_input_names.(lfuncs)), reduce(union, get_output_names.(lfuncs))
        #* Extract all parameters names of funcs and dfuncs
        ele_param_names = unique(reduce(union, get_param_names.(lfuncs)))
        #* Setup the name information of the hydrobucket
        infos = (input=ele_input_names, output=ele_output_names, param=ele_param_names)

        return new(
            name,
            infos,
            lfuncs,
        )
    end
end

function (ele::LagBucket)(
    input::AbstractMatrix,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver
)
    #* Extract the initial state of the parameters and bucket in the pas variable
    params = pas[:params]
    lag_weights = [lfunc.lag_func(params[get_param_names(lfunc)[1]]) for lfunc in ele.lfuncs]

    function solve_lag_flux(input_vec, lag_weight)
        #* 首先将lagflux转换为discrete problem
        function lag_prob(u, p, t)
            u = circshift(u, -1)
            u[end] = 0.0
            input_vec[Int(t)] .* p[:weight] .+ u
        end

        prob = DiscreteProblem(lag_prob, lag_weight, (timeidx[1], timeidx[end]),
            ComponentVector(weight=lag_weight))
        #* 求解这个问题
        sol = solve(prob, FunctionMap())
        sol
    end

    sols = solve_lag_flux.(eachslice(input, dims=1), lag_weights)
    reduce(hcat, [Array(sol)[1, :] for sol in sols])'
end

function run_multi_fluxes(
    ele::LagBucket;
    input::AbstractArray,
    pas::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and bucket in the pas variable
    #* var_name * [weight_len * node_num]
    params = pas[:params]
    lag_weights = [[lfunc.lag_func(params[ptype][get_param_names(lfunc)[1]]) for ptype in ptypes] for lfunc in ele.lfuncs]

    function solve_lag_flux(input_vec, lag_weight)
        #* 首先将lagflux转换为discrete problem
        function lag_prob(u, p, t)
            tmp_u = circshift(u, -1)
            tmp_u[end] = 0.0
            input_vec[Int(t)] .* p[:w] .+ tmp_u
        end
        prob = DiscreteProblem(lag_prob, lag_weight, (1, length(input_vec)), ComponentVector(w=lag_weight))
        sol = solve(prob, FunctionMap())
        Array(sol)[1, :]
    end

    sols = map(eachindex(ele.lfuncs)) do idx
        node_sols = reduce(hcat, solve_lag_flux.(eachslice(input[idx, :, :], dims=1), lag_weights[idx]))
        node_sols
    end
    if length(sols) > 1
        sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
        return permutedims(sol_arr, (3, 1, 2))
    else
        return reshape(sols[1], 1, size(input)[3], size(input)[2])
    end
end
