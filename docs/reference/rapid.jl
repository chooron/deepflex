
"""
    VectorRoute <: AbstractVectorRoute

A structure representing a vector-based routing scheme for hydrological modeling.

# Fields
- `rfunc::AbstractVector{<:AbstractRouteFlux}`: Vector of routing flux functions for each node.
- `network::DiGraph`: A directed graph representing the routing network topology.
- `infos::NamedTuple`: Contains information about the VectorRoute instance, including input, output, state, and parameter names.

# Constructor
    VectorRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        network::DiGraph
    )

Constructs a `VectorRoute` object with the given name, routing flux function, and network structure.

# Arguments
- `name::Symbol`: A symbol representing the name of the routing scheme.
- `rfunc::AbstractRouteFlux`: The routing flux function to be applied at each node.
- `network::DiGraph`: A directed graph representing the routing network topology.

The constructor extracts variable names, parameter names, and neural network names from the provided
routing flux function to set up the internal information structure of the `VectorRoute` object.

Note: 来源于Rapid汇流模型
"""
struct RapidRoute <: AbstractRapidRoute
    "Routing adjacency matrix"
    adjacency::AbstractMatrix
    "grid subarea information, km2"
    subareas::AbstractVector
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function RapidRoute(;
        network::DiGraph,
        subareas::Union{AbstractVector,Number},
        name::Union{Symbol,Nothing}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_names)), :_vector_route) : name
        meta = HydroMeta(name=route_name, inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* generate adjacency matrix from network
        adjacency = adjacency_matrix(network)'
        #* Convert subareas to a vector if it's a single value
        subareas = subareas isa AbstractVector ? subareas : fill(subareas, nv(network))
        @assert length(subareas) == nv(network) "The length of subareas must be the same as the number of nodes, but got subareas: $(length(subareas)) and nodes: $(nv(network))"
        return new(
            adjacency,
            subareas,
            meta,
        )
    end
end

function (route::RapidRoute)(
    input::Array,
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    interp = get(config, :interp, LinearInterpolation)
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))
    delta_t = get(config, :delta_t, 1.0)
    solver = DiscreteSolver(alg=FunctionMap{true}())

    @assert all(ptype in keys(pas[:params]) for ptype in ptypes) "Missing required parameters. Expected all of $(keys(pas[:params])), but got $(ptypes)."

    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = interp.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    #* 计算出每个节点的面积转换系数
    area_coeffs = @. 24 * 3600 / (route.subareas * 1e6) * 1e3

    #* prepare the parameters for the routing function
    #* 参数可能存在转换需求,其中包括A通常是固定值
    k_ps = [pas[:params][ptype][:k] for ptype in ptypes]
    x_ps = [pas[:params][ptype][:x] for ptype in ptypes]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    A = (p) -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

    function route_ode!(du, u, p, t)
        q_gen = [itp_func(t) for itp_func in itp_funcs] .* area_coeffs
        #* Ax = b, x is the q_out(t+1)
        rflux_b = p.c0 * q_gen + p.c1 * (route.adjacency * q_out_t1 + q_gen) + p.c2 * q_out_t1
        #* solve the linear equation (simple solve by matrix inversion)
        du[:] = A(p) \ (rflux_b .- A(p) * u)
    end

    #* solve the ode
    sol_arr = solver(route_ode!, ComponentVector(c0=c0, c1=c1, c2=c2), zeros(size(input_mat)[1]), timeidx, convert_to_array=true)
    return reshape(sol_arr, 1, size(sol_arr)...)
end