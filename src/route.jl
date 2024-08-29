(flux::AbstractRouteFlux)(input::Vector, pas::ComponentVector; kwargs...) = error("This struct does not support Vector input in $(typeof(flux)) subtype of the AbstractRouteFlux")
(flux::AbstractRouteFlux)(input::Matrix, pas::ComponentVector; kwargs...) = error("This struct does not support Matrix input in $(typeof(flux)) subtype of the AbstractRouteFlux")

"""
    WeightSumRoute <: AbstractRouteFlux

Represents a weighted cumulative sum routing structure for hydrological modeling.

# Fields
- `infos::NamedTuple`: A named tuple containing routing information with keys:
  - `name::Symbol`: The name of the routing component.
  - `input::Symbol`: The symbol representing the input variable.
  - `output::Symbol`: The symbol representing the output variable.
  - `param::Symbol`: The symbol representing the weight parameter.

# Constructor
    WeightSumRoute(
        name::Symbol;
        input::Num,
        output::Num,
        param::Num
    )

Constructs a WeightSumRoute instance.

# Arguments
- `name::Symbol`: The name of the routing component.
- `input::Num`: The input variable.
- `output::Num`: The output variable.
- `param::Num`: The weight parameter.

# Description
WeightSumRoute applies a weighted cumulative sum operation to the input,
where the weights are specified by the `param` parameter. This routing method
is useful for scenarios where the contribution of each input node needs to be
weighted differently in the cumulative output.
"""
struct WeightSumRoute <: AbstractRouteFlux
    "Routing information: keys contain: input, output, param, state, nn"
    infos::NamedTuple

    function WeightSumRoute(
        name::Symbol;
        input::Num,
        output::Num,
        param::Num,
    )
        input_name = Symbolic.tosymbol(input)
        output_name = Symbolic.tosymbol(output)
        param_name = Symbolic.tosymbol(param)
        infos = (name=name, input=input_name, output=output_name, param=param_name, state=Symbol[], nn=Symbol[])
        return new(
            infos,
        )
    end
end

function (route::WeightSumRoute)(
    input::AbstractArray,
    pas::ComponentVector;
    ptypes::AbstractVector{Symbol},
    kwargs...
)
    input_mat = input[1, :, :]
    weight_params = [pas[:params][ptype][route.infos[:param]] for ptype in ptypes]
    weight_result = sum(input_mat .* weight_params, dims=1)
    # expand dims
    output_arr = reduce(vcat, repeat(weight_result, size(input_mat)[1]))
    reshape(1, size(output_arr)...)
end


"""
    GridRoute(name::Symbol; rfunc::AbstractRouteFlux, flwdir::AbstractMatrix, positions::AbstractVector)

Represents a grid-based routing structure for hydrological modeling.

# Arguments
- `name::Symbol`: A symbol representing the name of the GridRoute instance.
- `rfunc::AbstractRouteFlux`: The routing function used for flow calculations.
- `flwdir::AbstractMatrix`: A matrix representing the flow direction for each grid cell.
- `positions::AbstractVector`: A vector of positions for each node in the grid.

# Fields
- `rfunc::AbstractRouteFlux`: The routing function used for flow calculations.
- `flwdir::AbstractMatrix`: A matrix representing the flow direction for each grid cell.
- `positions::AbstractVector`: A vector of positions for each node in the grid.
- `infos::NamedTuple`: Contains information about the GridRoute instance, including input, output, state, and parameter names.

# Description
GridRoute is a structure that represents a grid-based routing system in a hydrological model. 
It uses a specified routing function (`rfunc`) to calculate flow between grid cells based on 
the provided flow direction matrix (`flwdir`) and node positions (`positions`).

The `infos` field stores metadata about the GridRoute instance, including names of inputs, 
outputs, states, parameters, and neural networks (if applicable) derived from the routing function.

This structure is particularly useful for modeling water flow across a landscape represented as a grid, 
where each cell has a specific flow direction and contributes to downstream flow.
"""
struct GridRoute <: AbstractGridRoute
    "Routing function"
    rfunc::AbstractRouteFlux
    "Flow direction matrix"
    flwdir::AbstractMatrix
    "Node position information"
    positions::AbstractVector
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function GridRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        flwdir::AbstractMatrix,
        positions::AbstractVector,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            rfunc,
            flwdir,
            positions,
            infos,
        )
    end
end

"""
input dims: node_num * ts_len
"""
function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
    #* 转换为input的稀疏矩阵
    input_arr = Array(sparse(
        [pos[1] for pos in positions],
        [pos[2] for pos in positions],
        input,
        size(flwdir)[1], size(flwdir)[2]
    ))
    #* 计算权重求和结果
    input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg)
                                for (code, arg) in zip(d8_codes, d8_nn_pads)]))
    #* 裁剪输入矩阵边框
    clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
    #* 将输入矩阵转换为向量
    collect([clip_arr[pos[1], pos[2]] for pos in positions])
end

function (route::GridRoute)(
    input::AbstractArray,
    pas::ComponentVector;
    timeidx::AbstractVector,
    ptypes::AbstractVector{Symbol},
    kwargs...
)
    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = LinearInterpolation.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ptypes)
    flux_initstates = get_rflux_initstates(route.rfunc; pas, ptypes)

    function grid_route_ode!(du, u, p, t)
        q_in = u[:q_in]
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], q_gen, p)
        new_q_in = grid_routing(q_out, route.positions, route.flwdir)
        du[:q_in] = new_q_in .- q_in
    end

    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    prob = ODEProblem(grid_route_ode!, init_states, (1, size(input_mat)[2]), pas[:params])
    sol = solve(prob, Tsit5(), saveat=timeidx)

    # Extract flux_states and q_in for each time step
    flux_states_matrix = reduce(hcat, [u.flux_states for u in sol.u])
    q_in_matrix = reduce(hcat, [u.q_in for u in sol.u])

    q_out = cal_flux_q_out.(eachslice(flux_states_matrix, dims=2), eachslice(q_in_matrix, dims=2), eachslice(input_mat, dims=2), Ref(pas[:params]))
    q_out_mat = reduce(hcat, q_out)
    # Convert q_out_mat to 1 x mat size
    q_out_reshaped = reshape(q_out_mat, 1, size(q_out_mat)...)
    return q_out_reshaped
end

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
"""
struct VectorRoute <: AbstractVectorRoute
    "Routing function"
    rfunc::AbstractRouteFlux
    "Routing network"
    network::DiGraph
    "Routing information: keys contain: input, output, param, state, nn"
    infos::NamedTuple

    function VectorRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        network::DiGraph,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(
            rfunc,
            network,
            infos,
        )
    end
end

"""
step route
"""
function (route::VectorRoute)(
    input::AbstractArray,
    pas::ComponentVector;
    timeidx::AbstractVector,
    ptypes::AbstractVector{Symbol},
    kwargs...
)
    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = LinearInterpolation.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ptypes)
    flux_initstates = get_rflux_initstates(route.rfunc; pas, ptypes)

    sorted_idxes = topological_sort(route.network)
    up_idxes = [inneighbors(route.network, cur_idx) for cur_idx in sorted_idxes]

    function vec_route_ode!(du, u, p, t)
        q_in = u[:q_in]
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], q_gen, p)
        for (sorted_idx, up_idx) in zip(sorted_idxes, up_idxes)
            q_out[sorted_idx] = sum(q_out[up_idx]) .+ q_out[sorted_idx]
        end
        du[:q_in] = q_out .- q_in
    end

    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    prob = DiscreteProblem(vec_route_ode!, init_states, (1, size(input_mat)[2]), pas[:params])
    sol = solve(prob, FunctionMap())

    flux_states_matrix = reduce(hcat, [u.flux_states for u in sol.u])
    q_in_matrix = reduce(hcat, [u.q_in for u in sol.u])

    q_out = cal_flux_q_out.(eachslice(flux_states_matrix, dims=2), eachslice(q_in_matrix, dims=2), eachslice(input_mat, dims=2), Ref(pas[:params]))
    q_out
end