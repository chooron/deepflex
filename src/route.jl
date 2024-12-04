"""
    HydroRoute(; rfunc::AbstractHydroFlux, rstate::Num, projfunc::AbstractHydroFlux, name::Union{Symbol,Nothing}=nothing)

Represents a routing structure for hydrological modeling.

# Arguments
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `rstate::Num`: The state variable for routing.
- `projfunc::AbstractHydroFlux`: Function for projecting outflow to downstream nodes.
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance. If not provided, will be automatically generated.

# Fields
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `projfunc::Function`: Function for projecting outflow to downstream nodes.
- `meta::HydroMeta`: Contains metadata about the routing instance, including input, output, state, parameter and neural network names.

# Description
HydroRoute is a structure that represents a routing system in a hydrological model.
It uses a specified routing function (`rfunc`) to calculate flow between nodes and a 
projection function (`projfunc`) to determine how water moves between connected nodes.

The routing process involves:
1. Calculating outflow from each node using the routing function
2. Projecting outflows to downstream nodes using the projection function
3. Updating node states based on inflow, outflow and locally generated runoff

The structure supports both traditional parameter-based routing functions and neural network
based routing functions through the AbstractHydroFlux interface.

The metadata (`meta`) is automatically constructed from the provided functions and contains:
- Input names (excluding the routing state variable)
- Output names
- Parameter names
- State name (from `rstate`)
- Neural network names (if any)

This structure serves as the base for more specific routing implementations like GridRoute
and VectorRoute.
"""
struct HydroRoute{F<:AbstractHydroFlux,PF<:Function,M<:HydroMeta} <: AbstractHydroRoute
    "Routing function"
    rfunc::F
    "Outflow projection function"
    projfunc::PF
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::M

    function HydroRoute(;
        rfunc::F,
        rstate::Num,
        projfunc::PF,
        name::Union{Symbol,Nothing}=nothing,
    ) where {F<:AbstractHydroFlux,PF<:Function}
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names = get_var_names(rfunc)
        state_name = Symbolics.tosymbol(rstate)
        input_names = setdiff(input_names, [state_name])
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_names)), :_grid_route) : name
        meta = HydroMeta(route_name, input_names, output_names, param_names, [state_name], nn_names)

        return new{F,PF,typeof(meta)}(
            rfunc,
            projfunc,
            meta,
        )
    end
end

"""
    GridRoute(;
        rfunc::AbstractHydroFlux,
        rstate::Num,
        flwdir::AbstractMatrix,
        positions::AbstractVector,
        aggtype::Symbol=:type1,
        name::Union{Symbol,Nothing}=nothing
    )

Create a HydroRoute instance for grid-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `flwdir::AbstractMatrix`: Flow direction matrix using D8 encoding (1-128)
- `positions::AbstractVector`: Vector of (row,col) positions for each node in the grid
- `aggtype::Symbol=:matrix`: Aggregation type for flow routing:
    - `:matrix`: Uses graph-based adjacency matrix
    - `:network`: Uses network-based adjacency matrix
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance

# Returns
`HydroRoute`: A configured routing instance for grid-based networks

# Description
Creates a routing structure for grid-based river networks using either a matrix-based
approach (matrix) or network-based approach (network) for flow accumulation. The function
verifies that node positions match node IDs and constructs appropriate projection
functions based on the chosen aggregation type.
"""
function GridRoute(;
    rfunc::AbstractHydroFlux,
    rstate::Num,
    flwdir::AbstractMatrix,
    positions::AbstractVector,
    aggtype::Symbol=:matrix,
    name::Union{Symbol,Nothing}=nothing
)
    if aggtype == :matrix
        d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
        d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]

        #* input dims: node_num * ts_len
        function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
            #* 转换为input的稀疏矩阵
            input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
            #* 计算权重求和结果
            input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
            #* 裁剪输入矩阵边框
            clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
            #* 将输入矩阵转换为向量
            collect([clip_arr[pos[1], pos[2]] for pos in positions])
        end
        #* build the outflow projection function
        projfunc = (outflow) -> grid_routing(outflow, positions, flwdir)

    elseif aggtype == :network
        network = build_grid_digraph(flwdir, positions)
        #* build the outflow projection function
        adjacency = adjacency_matrix(network)'
        projfunc = (outflow) -> adjacency * outflow
    else
        @error "the $aggtype is not support"
    end

    return HydroRoute(; rfunc, rstate, projfunc, name)
end

"""
    VectorRoute(;
        rfunc::AbstractHydroFlux,
        rstate::Num,
        network::DiGraph,
        name::Union{Symbol,Nothing}=nothing
    )

Create a HydroRoute instance for vector-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `network::DiGraph`: Directed graph representing the river network connectivity
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance

# Returns
`HydroRoute`: A configured routing instance for vector-based networks

# Description
Creates a routing structure for vector-based river networks using a graph-based approach
for flow accumulation. The function verifies that the number of nodes matches the number
of node IDs and constructs a projection function based on the network's adjacency matrix.
The adjacency matrix is used to route flow between connected nodes in the network.
"""
function VectorRoute(;
    rfunc::AbstractHydroFlux,
    rstate::Num,
    network::DiGraph,
    name::Union{Symbol,Nothing}=nothing,
)
    #* generate adjacency matrix from network
    adjacency = adjacency_matrix(network)'
    #* build the outflow projection function
    projfunc = (outflow) -> adjacency * outflow
    return HydroRoute(; rfunc, rstate, projfunc, name)
end

function _get_parameter_extractors(route::HydroRoute, pas::ComponentVector, ptypes::AbstractVector{Symbol})
    if route.rfunc isa AbstractNeuralFlux
        @assert all(nn_name in keys(pas[:nn]) for nn_name in get_nn_names(route)) "Missing required neural networks. Expected all of $(get_nn_names(ele)), but got $(keys(pas[:nn]))."
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in get_nn_names(route.rfunc)]
        param_func = (p) -> Ref([p[:nn][idx] for idx in nn_params_idx])
    else
        for ptype in ptypes
            @assert all(param_name in keys(pas[:params][ptype]) for param_name in get_param_names(route.rfunc)) "Missing required parameters. Expected all of $(get_param_names(ele)), but got $(keys(pas[:params][ptype])) at param type: $ptype."
        end
        rflux_params_idx = [getaxes(pas[:params][ptypes[1]])[1][nm].idx for nm in get_param_names(route.rfunc)]
        param_func = (p) -> [p[:params][ptype][rflux_params_idx] for ptype in ptypes]
    end
    return param_func
end

"""
    (route::HydroRoute)(input::Array, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Run the routing model for given input and parameters.

# Arguments
- `input`: Input data array with dimensions (variables × nodes × time)
- `pas`: ComponentVector containing parameters, initial states and neural network parameters (if applicable)
- `config`: Configuration options including:
  - `ptypes`: Parameter types to use for each node (default: all parameter types)
  - `interp`: Interpolation method for input data (default: LinearInterpolation)
  - `solver`: AbstractHydroSolver to use for ODE solving (default: ODESolver())
  - `timeidx`: Vector of time indices (default: 1:size(input,3))

# Returns
Array with dimensions (states+outputs × nodes × time) containing:
- Routing states for each node over time
- Routed outflow for each node over time

# Description
This function executes the routing model by:
1. Validating input dimensions and parameter/state configurations
2. Setting up interpolation and solver configurations
3. Building parameter functions for either neural network or regular routing
4. Solving the routing equations using the specified solver
5. Computing outflows based on states and parameters

The routing can use either neural network based routing functions (AbstractNeuralFlux) or
regular routing functions, with parameters extracted accordingly.
"""
function (route::HydroRoute{F,PF,M})(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
) where {F<:AbstractHydroFlux,PF<:Function,M<:HydroMeta,T<:Number}
    #* get the parameter types and state types
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    stypes = get(config, :stypes, collect(keys(pas[:initstates])))
    #* get the interpolation type and solver type
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ODESolver())
    #* get the time index
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))
    #* get the solver
    solver = get(config, :solver, ODESolver())

    #* check the length of ptypes
    @assert(length(ptypes) == size(input, 2),
        "The length of ptypes must match the number of nodes ($(size(input,2))), but got $(length(ptypes)) ptypes")
    @assert(all(stype in keys(pas[:initstates]) for stype in stypes),
        "Invalid HRU names. All names must be one of $(keys(pas[:initstates])), but got $(stypes)")

    #* Extract the idx range of each variable in params, this extraction method is significantly more efficient than extracting by name
    param_func = _get_parameter_extractors(route, pas, ptypes)
    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itp_funcs = interp.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)
    #* prepare the initial states matrix (dims: state_num * node_num)
    init_states_mat = reduce(hcat, [collect(pas[:initstates][stype][get_state_names(route)]) for stype in stypes])'

    function du_func(u, p, t)
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        q_out_vec = route.rfunc.func.(eachslice(hcat(q_gen, u), dims=1), param_func(p), Ref(t))
        q_out = reduce(vcat, q_out_vec)
        q_in = route.projfunc(q_out)
        q_in .+ q_gen .- q_out
    end

    #* solve the problem
    sol_arr = solver(du_func, pas, init_states_mat, timeidx, convert_to_array=true)
    sol_arr_permuted = permutedims(sol_arr, (2, 1, 3))
    cont_arr = cat(input, sol_arr_permuted, dims=1)
    output_vec = [route.rfunc.func.(eachslice(cont_arr[:, :, i], dims=2), param_func(pas), timeidx[i]) for i in 1:size(input)[3]]
    out_arr = reduce(hcat, reduce.(vcat, output_vec))
    #* return route_states and q_out
    return cat(sol_arr_permuted, reshape(out_arr, 1, size(out_arr)...), dims=1)
end

function (route::HydroRoute)(input::Vector{<:NamedTuple}, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)
    for i in eachindex(input)
        @assert all(input_name in keys(input[i]) for input_name in get_input_names(route)) "Missing required inputs. Expected all of $(get_input_names(route)), but got $(keys(input[i])) at $i input."
    end
    input_mats = [reduce(hcat, collect(input[i][k] for k in get_input_names(route))) for i in eachindex(input)]
    input_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), input_mats)
    route(input_arr, pas; config=config, kwargs...)
end


"""
    RapidRoute<: AbstractRoute

A structure representing a vector-based routing scheme for hydrological modeling.

# Fields
- `rfunc::AbstractVector{<:AbstractRouteFlux}`: Vector of routing flux functions for each node.
- `network::DiGraph`: A directed graph representing the routing network topology.
- `infos::NamedTuple`: Contains information about the VectorRoute instance, including input, output, state, and parameter names.

# Constructor
    RapidRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        network::DiGraph
    )

Constructs a `RapidRoute` object with the given name, routing flux function, and network structure.

# Arguments
- `name::Symbol`: A symbol representing the name of the routing scheme.
- `rfunc::AbstractRouteFlux`: The routing flux function to be applied at each node.
- `network::DiGraph`: A directed graph representing the routing network topology.

The constructor extracts variable names, parameter names, and neural network names from the provided
routing flux function to set up the internal information structure of the `RapidRoute` object.

Note: from Rapid
"""
struct RapidRoute <: AbstractRapidRoute
    "Routing adjacency matrix"
    adjacency::AbstractMatrix
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function RapidRoute(
        fluxes::Pair{Vector{Num},Vector{Num}};
        network::DiGraph,
        name::Union{Symbol,Nothing}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        inputs, outputs = fluxes[1], fluxes[2]
        @assert length(inputs) == length(outputs) == 1 "The length of inputs and outputs must be the 1, but got inputs: $(length(inputs)) and outputs: $(length(outputs))"
        input_names = Symbolics.tosymbol.(inputs)
        output_names = Symbolics.tosymbol.(outputs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = [:rapid_k, :rapid_x]
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_names)), :_vector_route) : name
        meta = HydroMeta(route_name, input_names, output_names, param_names, Symbol[], Symbol[])
        #* generate adjacency matrix from network
        adjacency = adjacency_matrix(network)'
        return new(
            adjacency,
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
    solver = get(config, :solver, DiscreteSolver())

    @assert all(ptype in keys(pas[:params]) for ptype in ptypes) "Missing required parameters. Expected all of $(keys(pas[:params])), but got $(ptypes)."

    #* var num * node num * ts len
    itp_funcs = interp.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)

    #* prepare the parameters for the routing function
    k_ps = [pas[:params][ptype][:k] for ptype in ptypes]
    x_ps = [pas[:params][ptype][:x] for ptype in ptypes]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    A = (p) -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

    function du_func(u, p, t)
        q_out_t1 = u
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        #* Ax = b, x is the q_out(t+1)
        rflux_b = p.c0 .* q_gen .+ p.c1 .* (route.adjacency * q_out_t1 .+ q_gen) .+ p.c2 .* q_out_t1
        #* solve the linear equation (simple solve by matrix inversion)
        A(p) \ (rflux_b .- A(p) * u)
    end

    #* solve the ode
    sol_arr = solver(du_func, ComponentVector(c0=c0, c1=c1, c2=c2), zeros(size(input)[2]), timeidx, convert_to_array=true)
    return reshape(sol_arr, 1, size(sol_arr)...)
end