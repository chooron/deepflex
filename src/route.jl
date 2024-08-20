struct HydroRoute <: AbstractRoute
    "routement information: keys contains: input, output, param, state"
    nameinfo::NamedTuple
    "route types"
    routefunc::Function

    function HydroRoute(
        inputs::Vector{Num},
        routefunc::Function
    )
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = map(s -> Symbol(s, :_routed), input_names)
        #* Setup the name information of the hydroroutement
        nameinfo = (input=input_names, output=output_names, param=[:k, :x])

        return new(
            nameinfo,
            routefunc
        )
    end
end

function (route::HydroRoute)(
    input::AbstractMatrix,
    params::ComponentVector;
    timeidx::AbstractVector
)
    #* Extract the initial state of the parameters and routement in the pas variable
    sol_arrs = route.routefunc.(eachslice(input, dims=1), Ref(params), timeidx)
    reduce(hcat, sol_arrs)'
end

struct UnitHydroRoute <: AbstractRoute
    "routement information: keys contains: input, output, param, state"
    nameinfo::NamedTuple
    """
    Hydrological lag fluxes, 
    combined with ordinary hydrological flux to construct ordinary differential equations
    """
    uhfuncs::Vector

    function UnitHydroRoute(
        inputs::Vector{Num},
        params::Union{Num,Vector{Num}},
        uhfuncs::Union{Function,Vector{Function}},
    )
        if !(params isa Vector)
            params = repeat([params], length(inputs))
        end
        if !(uhfuncs isa Vector)
            uhfuncs = repeat([uhfuncs], length(inputs))
        end
        @assert length(params) == length(uhfuncs) "number of the params should be equal with the number of uhfuncs"
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = map(s -> Symbol(s, :_routed), input_names)
        param_names = Symbolics.tosymbol.(params, escape=false)
        #* Setup the name information of the hydroroutement
        nameinfo = (input=input_names, output=output_names, param=param_names)

        return new(
            nameinfo,
            uhfuncs,
        )
    end
end

function (route::UnitHydroRoute)(
    input::AbstractMatrix,
    pas::ComponentVector;
)
    #* Extract the initial state of the parameters and routement in the pas variable
    params = pas[:params]
    lag_weights = [uhfunc(params[pname]) for (uhfunc, pname) in zip(route.uhfuncs, route.nameinfo[:param])]
    sol_arrs = solve_uhfunc.(eachslice(input, dims=1), lag_weights)
    reduce(hcat, sol_arrs)'
end

function run_multi_fluxes(
    route::UnitHydroRoute;
    input::AbstractArray,
    params::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(params)),
)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    #* var_name * [weight_len * node_num]
    lag_weights = [[uhfunc(params[ptype][pname]) for ptype in ptypes] for (uhfunc, pname) in zip(route.uhfuncs, route.nameinfo[:param])]

    sols = map(eachindex(route.lfuncs)) do idx
        node_sols = reduce(hcat, solve_uhfunc.(eachslice(input[idx, :, :], dims=1), lag_weights[idx]))
        node_sols
    end
    if length(sols) > 1
        sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
        return permutedims(sol_arr, (3, 1, 2))
    else
        return reshape(sols[1], 1, size(input)[3], size(input)[2])
    end
end

function run_multi_fluxes(
    route::HydroRoute;
    input::AbstractArray,
    params::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(params)),
)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    #* var_name * [weight_len * node_num]
    pytype_params = [params[ptype] for ptype in ptypes]
    sols = map(eachindex(ptypes)) do (idx)
        node_sols = reduce(hcat, route.routefunc.(eachslice(input[idx, :, :], dims=1), pytype_params[idx]))
        node_sols
    end
    if length(sols) > 1
        sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
        return permutedims(sol_arr, (3, 1, 2))
    else
        return reshape(sols[1], 1, size(input)[3], size(input)[2])
    end
end


function solve_mskfunc(input_vec, params)
    k, x, dt = params.k, params.x, params.dt
    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))
    function msk_prob(u, p, t)
        println(t)
        q0 = u[1]
        c0, c1, c2 = p
        input1 = input_vec[Int(t)]
        input0 = input_vec[Int(t)-1]
        new_q = (c0 * input1) + (c1 * input0) + (c2 * q0)
        [new_q]
    end

    prob = DiscreteProblem(msk_prob, [input_vec[1]], (2, length(input_vec)), ComponentVector(c0=c0, c1=c1, c2=c2))
    sol = solve(prob, FunctionMap())
    reduce(vcat, sol.u)
end
