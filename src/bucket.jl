"""
    HydroBucket(name::Symbol; funcs::Vector, dfuncs::Vector=StateFlux[])

Represents a hydrological bucket model component.

# Arguments
- `name::Symbol`: A symbol representing the name of the HydroBucket instance.
- `funcs::Vector`: A vector of flux functions that describe the hydrological processes.
- `dfuncs::Vector`: A vector of state derivative functions (default is an empty vector of StateFlux).

# Fields
- `flux_func::Function`: Combined function for calculating fluxes.
- `ode_func::Union{Nothing,Function}`: Function for ordinary differential equations (ODE) calculations, or nothing if not applicable.
- `infos::NamedTuple`: Contains metadata about the bucket, including input, output, state, parameter, and neural network names.

# Description
HydroBucket is a structure that encapsulates the behavior of a hydrological bucket model. 
It combines multiple flux functions and state derivative functions to model water movement 
and storage within a hydrological unit.

The structure automatically extracts relevant information from the provided functions to 
populate the `infos` field, which includes names of inputs, outputs, states, parameters, 
and neural networks (if applicable).

The `flux_func` and `ode_func` are constructed based on the provided `funcs` and `dfuncs`, 
enabling efficient calculation of fluxes and state changes over time.

This structure is particularly useful for building complex hydrological models by combining 
multiple HydroBucket instances to represent different components of a water system.

"""
struct HydroBucket <: AbstractBucket
    """
    Vector of flux functions describing hydrological processes.
    """
    funcs::Vector{<:AbstractFlux}
    """
    Vector of state derivative functions for ODE calculations.
    """
    dfuncs::Vector{<:AbstractStateFlux}
    """
    Combined function for calculating all hydrological fluxes.
    """
    flux_func::Function
    """
    Function for ordinary differential equations (ODE) calculations.
    Can be `nothing` if no ODE calculations are needed.
    """
    ode_func::Union{Nothing,Function}
    """
    Metadata about the bucket, including:
    - name: Symbol representing the bucket's name
    - input: Vector of input variable names
    - output: Vector of output variable names
    - state: Vector of state variable names
    - param: Vector of parameter names
    - nn: Vector of neural network names (if applicable)
    """
    meta::HydroMeta

    function HydroBucket(;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
        name::Union{Symbol, Nothing}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(funcs, dfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(vcat(funcs, dfuncs))
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(funcs)
        #* Setup the name information of the hydrobucket
        bucket_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), state_names)), :_bucket) : name
        meta = HydroMeta(name=bucket_name, inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(funcs, dfuncs, meta)

        return new(
            funcs,
            dfuncs,
            flux_func,
            ode_func,
            meta,
        )
    end
end

"""
    (ele::HydroBucket)(input::Matrix, pas::ComponentVector; timeidx::Vector, solver::AbstractSolver=ODESolver())
    (ele::HydroBucket)(input::Array, pas::ComponentVector; timeidx::Vector, ptypes::AbstractVector{Symbol}, solver::AbstractSolver=ODESolver())

Run the HydroBucket model for given input and parameters.

# Arguments
- `input`: Input data. Can be a Matrix or an Array.
- `pas`: ComponentVector containing parameters and initial states.
- `timeidx`: Vector of time indices.
- `solver`: AbstractSolver to use for ODE solving. Defaults to ODESolver().
- `ptypes`: (Only for Array input) AbstractVector of Symbol specifying parameter types.

# Returns
A matrix (for single node, dims: var_names × ts_len) or array (for multiple nodes, dims: var_names × node_names × ts_len) containing the model outputs and (if applicable) solved states.

# Details
- For Matrix input: Assumes input dimensions are [variables, time].
- For Array input: Allows for more complex input structures, requires `ptypes`.
- If the bucket has an ODE function, it solves the states over time.
- Calculates fluxes using the model's flux function.
- Concatenates solved states (if any) with calculated fluxes for output.

"""

function (ele::HydroBucket)(
    input::Matrix,
    pas::ComponentVector,
    timeidx::Vector{<:Number}=collect(1:size(input, 2));
    config::NamedTuple=NamedTuple(),
    kwargs...
)
    #* get kwargs
    solver = get(config, :solver, ODESolver())
    interp = get(config, :interp, LinearInterpolation)
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)

    @assert size(input, 1) == length(get_input_names(ele)) "Input dimensions mismatch. Expected $(length(get_input_names(ele))) variables, got $(size(input, 1))."
    @assert size(input, 2) == length(timeidx) "Time steps mismatch. Expected $(length(timeidx)) time steps, got $(size(input, 2))."

    #* Extract the initial state of the parameters and bucket in the pas variable
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of bucket at the specified timeidx
        solved_states = solve_prob(ele, input, pas, timeidx, solver=solver, interp=interp)
        #* Store the solved bucket state in fluxes
        fluxes = cat(input, solved_states, dims=1)
    else
        fluxes = input
        solved_states = nothing
    end

    #* extract params and nn params
    # Check if all required parameter names are present in pas[:params]
    @assert all(param_name in keys(pas[:params]) for param_name in get_param_names(ele)) "Missing required parameters. Expected all of $(get_param_names(ele)), but got $(keys(pas[:params]))."

    # Check if all required neural network names are present in pas[:nn] (if any)
    if !isempty(get_nn_names(ele))
        @assert all(nn_name in keys(pas[:nn]) for nn_name in get_nn_names(ele)) "Missing required neural networks. Expected all of $(get_nn_names(ele)), but got $(keys(pas[:nn]))."
    end

    params_vec = collect([pas[:params][nm] for nm in get_param_names(ele)])
    nn_params_vec = !isempty(get_nn_names(ele)) ? collect(pas[:nn][nm] for nm in get_nn_names(ele)) : nothing
    #* calculate output, slice input on time dim, then calculate each output
    flux_output = ele.flux_func.(eachslice(fluxes, dims=2), Ref(params_vec), Ref(nn_params_vec), timeidx)
    #* convert vector{vector} to matrix
    flux_output_matrix = reduce(hcat, flux_output)

    #* merge output and state, if solved_states is not nothing, then cat it at the first dim
    output_matrix = isnothing(solved_states) ? flux_output_matrix : cat(solved_states, flux_output_matrix, dims=1)
    if convert_to_ntp
        return NamedTuple{Tuple(vcat(get_state_names(ele), get_output_names(ele)))}(eachslice(output_matrix, dims=1))
    else
        return output_matrix
    end
end

function (ele::HydroBucket)(
    input::Array,
    pas::ComponentVector,
    timeidx::Vector{<:Number}=collect(1:size(input, 3));
    config::NamedTuple=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)
    #* get kwargs
    solver = get(config, :solver, ODESolver())
    interp = get(config, :interp, LinearInterpolation)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)

    @assert size(input, 1) == length(get_input_names(ele)) "Input dimensions mismatch. Expected $(length(get_input_names(ele))) variables, got $(size(input, 1))."
    @assert size(input, 3) == length(timeidx) "Time steps mismatch. Expected $(length(timeidx)) time steps, got $(size(input, 3))."

    #* Extract the initial state of the parameters and bucket in the pas variable
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of bucket at the specified timeidx
        solved_states = solve_prob(ele, input, pas, timeidx, solver=solver, interp=interp)
        if solved_states == false
            solved_states = zeros(length(get_state_names(ele)), length(timeidx))
        end
        #* Store the solved bucket state in fluxes
        fluxes = cat(input, solved_states, dims=1)
    else
        fluxes = input
        solved_states = nothing
    end

    #* extract params and nn params
    params_collect = [pas[:params][ptype] for ptype in ptypes]
    #* check params input is correct
    for (ptype, params_item) in zip(ptypes, params_collect)
        @assert all(param_name in keys(params_item) for param_name in get_param_names(ele)) "Missing required parameters. Expected all of $(get_param_names(ele)), but got $(keys(params_item)) at param type: $ptype."
    end
    params_vec = collect([collect([params_item[pname] for pname in get_param_names(ele)]) for params_item in params_collect])

    # Check if all required neural network names are present in pas[:nn] (if any)
    if !isempty(get_nn_names(ele))
        @assert all(nn_name in keys(pas[:nn]) for nn_name in get_nn_names(ele)) "Missing required neural networks. Expected all of $(get_nn_names(ele)), but got $(keys(pas[:nn]))."
    end
    nn_params_vec = !isempty(get_nn_names(ele)) ? collect(pas[:nn][nm] for nm in get_nn_names(ele)) : nothing

    #* array dims: (num of node, sequence length, variable dim)
    ele_output_vec = [ele.flux_func.(eachslice(fluxes[:, i, :], dims=2), Ref(params_vec[i]), Ref(nn_params_vec), timeidx) for i in 1:size(input)[2]]
    ele_output_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, u) for u in ele_output_vec]), (1, 3, 2))

    #* merge state and output, if solved_states is not nothing, then cat it at the first dim
    final_output_arr = isnothing(solved_states) ? ele_output_arr : cat(solved_states, ele_output_arr, dims=1)
    final_output_arr
end

function (ele::HydroBucket)(input::NamedTuple, pas::ComponentVector, timeidx::Vector; config::NamedTuple=NamedTuple(), kwargs...)
    @assert all(input_name in keys(input) for input_name in get_input_names(ele)) "Missing required inputs. Expected all of $(get_input_names(ele)), but got $(keys(input))."
    for k in keys(input)
        @assert length(input[k]) == length(timeidx) "Time steps mismatch. Expected $(length(timeidx)) time steps, got $(length(input[k])) at input: $k."
    end
    input_matrix = Matrix(reduce(hcat, [input[k] for k in keys(input)])')
    ele(input_matrix, pas, timeidx; config=config, kwargs...)
end

function (ele::HydroBucket)(input::Vector{<:NamedTuple}, pas::ComponentVector, timeidx::Vector; config::NamedTuple=NamedTuple(), kwargs...)
    for i in eachindex(input)
        @assert all(input_name in keys(input[i]) for input_name in get_input_names(ele)) "Missing required inputs. Expected all of $(get_input_names(ele)), but got $(keys(input[i])) at $i input."
        @assert length(input[i]) == length(timeidx) "Time steps mismatch. Expected $(length(timeidx)) time steps, got $(length(input[i])) at $i input."
    end
    input_mats = [reduce(hcat, collect(input[i][k] for k in get_input_names(ele))) for i in eachindex(input)]
    input_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), input_mats)
    ele(input_arr, pas, timeidx; config=config, kwargs...)
end


"""
Solve the ordinary differential equations for a HydroBucket model.

This function handles two types of input arguments:

1. Single node input:
   - `input`: Matrix with dimensions (var_names × ts_len)
   - `pas`: ComponentVector with structure:
     ComponentVector(params=(p1=, p2=, ...), initstates=(...), nn=(...))

2. Multiple node input:
   - `input`: Array with dimensions (var_names × node_names × ts_len)
   - `pas`: ComponentVector with structure:
     ComponentVector(params=(node_1=(p1=, p2=, ...), node_2=(p1=, p2=, ...), ...), initstates=(...), nn=(...))

# Arguments
- `ele::HydroBucket`: The HydroBucket model instance
- `input`: Input data (Matrix for single node, dims: var_names × ts_len, Array for multiple nodes, dims: var_names × node_names × ts_len)
- `pas::ComponentVector`: Parameters and initial states
- `timeidx::Vector`: Time index vector
- `ptypes::Vector{Symbol}`: Parameter types (for multiple node input)
- `solver::AbstractSolver`: ODE solver to use (default: ODESolver())

# Returns
- `sol`: Solution of the ordinary differential equations, dims: state_names × ts_len (or state_names × node_names × ts_len for multiple nodes)
"""
function solve_prob(
    ele::HydroBucket,
    input::Matrix,
    pas::ComponentVector,
    timeidx::Vector{<:Number}=collect(1:size(input, 2));
    solver::AbstractSolver=ODESolver(),
    interp::Type{<:AbstractInterpolation}=LinearInterpolation,
)
    params, init_states = pas[:params], pas[:initstates]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_list = map(eachrow(input)) do var
        interp(var, timeidx, extrapolate=true)
    end
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params)[1][nm].idx for nm in get_param_names(ele)]
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    ode_param_func = (p) -> p[:params][params_idx]

    if !isempty(get_nn_names(ele))
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in get_nn_names(ele)]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function single_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        du[:] = ele.ode_func(ode_input, u, ode_params, nn_params, t)
    end

    #* Solve the problem using the solver wrapper
    sol = solver(single_ele_ode_func!, pas, collect(init_states[get_state_names(ele)]), timeidx)
    sol
end

function solve_prob(
    ele::HydroBucket,
    input::Array,
    pas::ComponentVector,
    timeidx::Vector{<:Number}=collect(1:size(input, 3));
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
    solver::AbstractSolver=ODESolver(),
    interp::Type{<:AbstractInterpolation}=LinearInterpolation,
)
    #* 针对多个相同的state function采用并行化计算,这样能够避免神经网络反复多次计算减少梯度计算反馈
    #* 同时将多组state function放到同一个ode function中,这种并行计算或能提高预测性能,
    #* 这样每个步长的输入维度就是:节点个数*输入变量数
    #* 当前只针对unit相同的同步求解:
    params, init_states = pas[:params], pas[:initstates]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_vecs = [interp.(eachslice(input[:, i, :], dims=1), Ref(timeidx), extrapolate=true) for i in 1:size(input)[2]]
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]

    #* 准备初始状态
    init_states_vec = collect([collect(init_states[ptype][get_state_names(ele)]) for ptype in ptypes])
    init_states_matrix = reduce(hcat, init_states_vec)

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params[ptypes[1]])[1][nm].idx for nm in get_param_names(ele)]
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    ode_param_func = (p) -> [p[:params][ptype][params_idx] for ptype in ptypes]

    #* 准备神经网络的参数
    if length(get_nn_names(ele)) > 0
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in get_nn_names(ele)]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function multi_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        tmp_output_vec = ele.ode_func.(ode_input, eachslice(u, dims=2), ode_params, Ref(nn_params), Ref(t))
        tmp_output = reduce(hcat, tmp_output_vec)
        du[:] = tmp_output
    end

    #* Solve the problem using the solver wrapper
    sol = solver(multi_ele_ode_func!, pas, init_states_matrix, timeidx)
    sol
end

