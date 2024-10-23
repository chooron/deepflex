(route::AbstractRoute)(input::Vector, pas::ComponentVector; kwargs::NamedTuple=NamedTuple()) = error("This struct does not support Vector input in $(typeof(route)) subtype of the AbstractRoute")
(route::AbstractRoute)(input::Matrix, pas::ComponentVector; kwargs::NamedTuple=NamedTuple()) = error("This struct does not support Matrix input in $(typeof(route)) subtype of the AbstractRoute")

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
struct WeightSumRoute <: AbstractSumRoute
    "routing function"
    rfunc::AbstractFlux
    "grid subarea information, km2"
    subareas::AbstractVector
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta
    # "output identifier" 
    #=     
    This part doesn't need to add an output id, because the matrix cannot be directly modified,
    so this part is all one value, representing the output result at the output id
    =#
    # outid::Symbol

    function WeightSumRoute(;
        rfunc::AbstractFlux=SimpleFlux([flow] => [flow_routed], exprs=[flow]),
        subareas::AbstractVector,
        name::Union{Symbol,Nothing}=nothing
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        push!(param_names, :route_weight)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), output_names)), :_weight_route) : name
        meta = HydroMeta(name=route_name, inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        return new(
            rfunc,
            subareas,
            meta,
        )
    end
end

function (route::WeightSumRoute)(
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::Vector{<:Number};
    kwargs...
)
    ptypes = get(kwargs, :ptypes, collect(keys(pas[:params])))
    #* 计算出每个节点的面积转换系数
    area_coefs = @. 24 * 3600 / (route.subareas * 1e6) * 1e3

    rfunc_output = route.rfunc(input, pas, timeidx, ptypes=ptypes)
    weight_params = [pas[:params][ptype][:route_weight] for ptype in ptypes]
    weight_result = sum(rfunc_output[1, :, :] .* weight_params .* area_coefs, dims=1)
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
    "grid subarea information, km2"
    subareas::AbstractVector
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function GridRoute(;
        rfunc::AbstractRouteFlux,
        flwdir::AbstractMatrix,
        positions::AbstractVector,
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
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_names)), :_grid_route) : name
        meta = HydroMeta(name=route_name, inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* Convert subareas to a vector if it's a single value
        subareas = subareas isa AbstractVector ? subareas : fill(subareas, length(positions))
        return new(
            rfunc,
            flwdir,
            positions,
            subareas,
            meta,
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
    pas::ComponentVector,
    timeidx::Vector{<:Number};
    config::NamedTuple=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ODESolver())

    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = interp.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    #* 计算出每个节点的面积转换系数
    area_coefs = @. 24 * 3600 / (route.subareas * 1e6) * 1e3

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ptypes)
    flux_initstates = get_rflux_initstates(route.rfunc; input=input_mat, pas=pas, ptypes=ptypes)

    function grid_route_ode!(du, u, p, t)
        # 提取单元产流
        q_gen = [itp_func(t) for itp_func in itp_funcs] .* area_coefs
        # 计算单元出流,更新单元出流状态
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], q_gen, p)
        # 计算出流的汇流结果
        new_q_in = grid_routing(q_out, route.positions, route.flwdir)
        # 更新状态
        du[:q_in] = new_q_in .- u[:q_in]
    end
    #* prepare init states
    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    #* solve the ode
    sol = solver(grid_route_ode!, pas[:params], init_states, timeidx, convert_to_array=false)
    if SciMLBase.successful_retcode(sol)
        # Extract flux_states and q_in for each time step
        flux_states_matrix = reduce(hcat, [u.flux_states for u in sol.u])
        q_in_matrix = reduce(hcat, [u.q_in for u in sol.u])
    else
        @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
        flux_states_matrix = zeros(size(flux_initstates)..., length(timeidx))
        q_in_matrix = zeros(size(input_mat)[1], length(timeidx))
    end

    q_out = cal_flux_q_out.(eachslice(flux_states_matrix, dims=2), eachslice(q_in_matrix, dims=2), eachslice(input_mat, dims=2), Ref(pas[:params]))
    q_out_mat = reduce(hcat, q_out)
    # Convert q_out_mat and flux_states_matrix to 1 x mat size
    q_out_reshaped = reshape(q_out_mat, 1, size(q_out_mat)...)
    flux_states_matrix_reshaped = reshape(flux_states_matrix, 1, size(flux_states_matrix)...)
    #  return flux_states and q_out
    return cat(flux_states_matrix_reshaped, q_out_reshaped, dims=1)
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
    "Routing adjacency matrix"
    adjacency::AbstractMatrix
    "grid subarea information, km2"
    subareas::AbstractVector
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function VectorRoute(;
        rfunc::AbstractRouteFlux,
        network::DiGraph,
        subareas::Union{AbstractVector,Number},
        name::Union{Symbol,Nothing}=nothing
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
        subareas = subareas isa AbstractVector ? subareas : fill(subareas, length(nv(network)))
        return new(
            rfunc,
            network,
            adjacency,
            subareas,
            meta,
        )
    end
end

"""
step route
"""
function (route::VectorRoute)(
    input::AbstractArray,
    pas::ComponentVector,
    timeidx::AbstractVector;
    config::NamedTuple=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ODESolver())

    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    input_mat = input[1, :, :]
    itp_funcs = interp.(eachslice(input_mat, dims=1), Ref(timeidx), extrapolate=true)

    #* 计算出每个节点的面积转换系数
    area_coefs = @. 24 * 3600 / (route.subareas * 1e6) * 1e3

    cal_flux_q_out!, cal_flux_q_out = get_rflux_func(route.rfunc; pas, ptypes)
    flux_initstates = get_rflux_initstates(route.rfunc; input=input_mat, pas=pas, ptypes=ptypes)

    function vec_route_ode!(du, u, p, t)
        q_in = u[:q_in]
        q_gen = [itp_func(t) for itp_func in itp_funcs] .* area_coefs
        q_out = cal_flux_q_out!(du, u[:flux_states], u[:q_in], q_gen, p)
        # update up_in
        q_in_updated = route.adjacency * q_out
        du[:q_in] = q_in_updated .- q_in
    end

    #* prepare init states
    init_states = ComponentVector(flux_states=flux_initstates, q_in=zeros(size(input_mat)[1]))
    #* solve the ode
    sol = solver(vec_route_ode!, pas[:params], init_states, timeidx, convert_to_array=false)
    if SciMLBase.successful_retcode(sol)
        # Extract flux_states and q_in for each time step
        flux_states_matrix = reduce(hcat, [u.flux_states for u in sol.u])
        q_in_matrix = reduce(hcat, [u.q_in for u in sol.u])
    else
        @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
        flux_states_matrix = zeros(size(flux_initstates)..., length(timeidx))
        q_in_matrix = zeros(size(input_mat)[1], length(timeidx))
    end

    q_out = cal_flux_q_out.(eachslice(flux_states_matrix, dims=2), eachslice(q_in_matrix, dims=2), eachslice(input_mat, dims=2), Ref(pas[:params]))
    q_out_mat = reduce(hcat, q_out)
    # Convert q_out_mat and flux_states_matrix to 1 x mat size
    q_out_reshaped = reshape(q_out_mat, 1, size(q_out_mat)...)
    flux_states_matrix_reshaped = reshape(flux_states_matrix, 1, size(flux_states_matrix)...)
    return cat(flux_states_matrix_reshaped, q_out_reshaped, dims=1)
end