"""
    HydroRoute(; rfunc::AbstractHydroFlux, rstate::Num, projfunc::AbstractHydroFlux, hrunames::Vector{Symbol}, name::Union{Symbol,Nothing}=nothing)

Represents a routing structure for hydrological modeling.

# Arguments
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `rstate::Num`: The state variable for routing.
- `projfunc::AbstractHydroFlux`: Function for projecting outflow to downstream nodes.
- `hrunames::Vector{Symbol}`: A vector of node identifiers.
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance. If not provided, will be automatically generated.

# Fields
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `projfunc::Function`: Function for projecting outflow to downstream nodes.
- `hrunames::Vector{Symbol}`: Node identifiers.
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
    "hru names"
    hrunames::Vector{Symbol}
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::M

    function HydroRoute(;
        rfunc::F,
        rstate::Num,
        projfunc::PF,
        hrunames::Vector{Symbol},
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
            hrunames,
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
        hrunames::Vector{Symbol},
        aggtype::Symbol=:type1,
        name::Union{Symbol,Nothing}=nothing
    )

Create a HydroRoute instance for grid-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `flwdir::AbstractMatrix`: Flow direction matrix using D8 encoding (1-128)
- `positions::AbstractVector`: Vector of (row,col) positions for each node in the grid
- `hrunames::Vector{Symbol}`: Vector of node identifiers
- `aggtype::Symbol=:type1`: Aggregation type for flow routing:
    - `:type1`: Uses D8 flow direction codes with padding
    - `:type2`: Uses graph-based adjacency matrix
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance

# Returns
`HydroRoute`: A configured routing instance for grid-based networks

# Description
Creates a routing structure for grid-based river networks using either a padding-based
approach (type1) or graph-based approach (type2) for flow accumulation. The function
verifies that node positions match node IDs and constructs appropriate projection
functions based on the chosen aggregation type.
"""
function GridRoute(;
    rfunc::AbstractHydroFlux,
    rstate::Num,
    flwdir::AbstractMatrix,
    positions::AbstractVector,
    hrunames::Vector{Symbol},
    aggtype::Symbol=:type1,
    name::Union{Symbol,Nothing}=nothing
)
    @assert length(positions) == length(hrunames) "The length of positions must be the same as the length of hrunames, but got positions: $(length(positions)) and hrunames: $(length(hrunames))"
    if aggtype == :type1
        d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
        d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]

        """
        input dims: node_num * ts_len
        """
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

    elseif aggtype == :type2
        network = build_grid_digraph(flwdir, positions)
        #* build the outflow projection function
        adjacency = adjacency_matrix(network)'
        projfunc = (outflow) -> adjacency * outflow
    else
        @error "the $aggtype is not support"
    end

    return HydroRoute(; rfunc, rstate, projfunc, hrunames, name)
end

"""
    VectorRoute(;
        rfunc::AbstractHydroFlux,
        rstate::Num,
        network::DiGraph,
        hrunames::Vector{Symbol},
        name::Union{Symbol,Nothing}=nothing
    )

Create a HydroRoute instance for vector-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `network::DiGraph`: Directed graph representing the river network connectivity
- `hrunames::Vector{Symbol}`: Vector of node identifiers
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
    hrunames::Vector{Symbol},
    name::Union{Symbol,Nothing}=nothing,
)
    @assert length(hrunames) == nv(network) "The length of hrunames must be the same as the number of nodes, but got hrunames: $(length(hrunames)) and nodes: $(nv(network))"
    #* generate adjacency matrix from network
    adjacency = adjacency_matrix(network)'
    #* build the outflow projection function
    projfunc = (outflow) -> adjacency * outflow
    return HydroRoute(; rfunc, rstate, projfunc, hrunames, name)
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
    stypes = get(config, :stypes, route.hrunames)
    #* get the interpolation type and solver type
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ODESolver())
    #* get the time index
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))

    #* check the length of ptypes, route.hrunames
    @assert(length(ptypes) == size(input, 2),
        "The length of ptypes must match the number of nodes ($(size(input,2))), but got $(length(ptypes)) ptypes")
    @assert(all(ptype in keys(pas[:params]) for ptype in ptypes),
        "Invalid parameter types. All ptypes must be one of $(keys(pas[:params])), but got $(ptypes)")
    @assert(all(stype in keys(pas[:initstates]) for stype in stypes),
        "Invalid HRU names. All names must be one of $(keys(pas[:initstates])), but got $(stypes)")
    @assert(all(stype in keys(pas[:initstates]) for stype in route.hrunames),
        "Invalid HRU names. All names must be one of $(keys(pas[:initstates])), but got $(route.hrunames)")   

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    if route.rfunc isa AbstractNeuralFlux
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in get_nn_names(route.rfunc)]
        param_func = (p) -> Ref([p[:nn][idx] for idx in nn_params_idx])
    else
        rflux_params_idx = [getaxes(pas[:params][ptypes[1]])[1][nm].idx for nm in get_param_names(route.rfunc)]
        param_func = (p) -> [p[:params][ptype][rflux_params_idx] for ptype in ptypes]
    end

    #* solve the problem
    sol_arr = solve_prob(route, input, pas, paramfunc=param_func, timeidx=timeidx, solver=solver, interp=interp)
    sol_arr_permuted = permutedims(sol_arr, (2, 1, 3))
    cont_arr = cat(input, sol_arr_permuted, dims=1)
    output_vec = [route.rfunc.func.(eachslice(cont_arr[:, :, i], dims=2), param_func(pas), timeidx[i]) for i in 1:size(input)[3]]
    out_arr = reduce(hcat, reduce.(vcat, output_vec))
    #  return route_states and q_out
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
    solve_prob(
        route::HydroRoute,
        input::Array,
        pas::ComponentVector;
        paramfunc::Function,
        timeidx::Vector{<:Number}=collect(1:size(input, 3)),
        solver::AbstractHydroSolver=ODESolver(),
        interp::Type{<:AbstractInterpolation}=LinearInterpolation,
    )

Solve the routing problem for a HydroRoute model.

# Arguments
- `route::HydroRoute`: The routing model
- `input::Array`: Input data array (dims: variables × nodes × time)
- `pas::ComponentVector`: Parameters and initial states
- `paramfunc::Function`: Function to extract parameters for the routing function
- `timeidx::Vector{<:Number}`: Vector of time indices (default: 1:size(input, 3))
- `solver::AbstractHydroSolver`: Solver to use for ODE solving (default: ODESolver())
- `interp::Type{<:AbstractInterpolation}`: Interpolation method for input data (default: LinearInterpolation)

# Returns
Array containing the solved states over time (dims: nodes × states × time)

# Description
This function solves the routing equations for a hydrological network. It:
1. Interpolates the input data for continuous time calculations
2. Prepares initial states for all nodes
3. Constructs and solves the coupled ODEs that represent routing between nodes
4. Returns the evolution of states over time

The routing equations consider both local runoff generation (from input) and flow between nodes
based on the routing function and projection matrix specified in the HydroRoute model.
"""

function solve_prob(
    route::HydroRoute{F,PF,M},
    input::AbstractArray{T,3},
    pas::ComponentVector;
    paramfunc::Function,
    timeidx::Vector{<:Number}=collect(1:size(input, 3)),
    solver::AbstractHydroSolver=ODESolver(),
    interp::Type{<:AbstractInterpolation}=LinearInterpolation,
) where {F<:AbstractHydroFlux,PF<:Function,M<:HydroMeta,T<:Number}
    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itp_funcs = interp.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)
    #* prepare the initial states matrix (dims: state_num * node_num)
    init_states_mat = reduce(hcat, [collect(pas[:initstates][hname][get_state_names(route)]) for hname in route.hrunames])'

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function route_ode!(du, u, p, t)
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        q_out_vec = route.rfunc.func.(eachslice(hcat(q_gen, u), dims=1), paramfunc(p), Ref(t))
        q_out = reduce(vcat, q_out_vec)
        q_in = route.projfunc(q_out)
        du[:] = q_in .+ q_gen .- q_out
    end

    #* Solve the problem using the solver wrapper
    sol = solver(route_ode!, pas, init_states_mat, timeidx, convert_to_array=true)
    sol
end