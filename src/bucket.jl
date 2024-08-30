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
    Hydrological flux functions
    """
    flux_func::Function
    """
    Hydrological ode functions
    """
    ode_func::Union{Nothing,Function}
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function HydroBucket(
        name::Symbol;
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(funcs, dfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(vcat(funcs, dfuncs))
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(funcs)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(funcs, dfuncs, infos)

        return new(
            flux_func,
            ode_func,
            infos,
        )
    end
end

function (ele::HydroBucket)(
    input::Matrix,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver(),
)
    #* Extract the initial state of the parameters and bucket in the pas variable
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of bucket at the specified timeidx
        solved_states = solve_single_prob(ele, input=input, pas=pas, timeidx=timeidx, solver=solver)
        if solved_states == false
            solved_states = zeros(length(ele.infos[:state]), length(timeidx))
        end
        #* Store the solved bucket state in fluxes
        fluxes = cat(input, solved_states, dims=1)
    else
        fluxes = input
        solved_states = nothing
    end

    #* excute other fluxes formula
    params_vec = collect([pas[:params][nm] for nm in ele.infos[:param]])
    if !isempty(ele.infos[:nn])
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.infos[:nn]])
    else
        nn_params_vec = nothing
    end
    flux_output = ele.flux_func.(eachslice(fluxes, dims=2), Ref(params_vec), Ref(nn_params_vec))
    #* convert vector{vector} to matrix
    flux_output_matrix = reduce(hcat, flux_output)

    #* merge output and state
    if isnothing(solved_states)
        output_matrix = flux_output_matrix
    else
        output_matrix = cat(solved_states, flux_output_matrix, dims=1)
    end
    output_matrix
end

function (ele::HydroBucket)(
    input::Array,
    pas::ComponentVector;
    timeidx::Vector,
    ptypes::AbstractVector{Symbol},
    solver::AbstractSolver=ODESolver(),
)
    #* Extract the initial state of the parameters and bucket in the pas variable
    if !isnothing(ele.ode_func)
        #* Call the solve_prob method to solve the state of bucket at the specified timeidx
        solved_states = solve_multi_prob(ele, input=input, pas=pas, timeidx=timeidx, solver=solver)
        if solved_states == false
            solved_states = zeros(length(ele.infos[:state]), length(timeidx))
        end
        #* Store the solved bucket state in fluxes
        fluxes = cat(input, solved_states, dims=1)
    else
        fluxes = input
        solved_states = nothing
    end

    params_vec = collect([collect([pas[:params][ptype][pname] for pname in ele.infos[:param]]) for ptype in ptypes])
    if !isempty(ele.infos[:nn])
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.infos[:nn]])
    else
        nn_params_vec = nothing
    end
    #* array dims: (num of node, sequence length, variable dim)
    ele_output_vec = [ele.flux_func.(eachslice(fluxes[:, i, :], dims=2), Ref(params_vec[i]), Ref(nn_params_vec)) for i in 1:size(input)[2]]
    ele_output_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, u) for u in ele_output_vec]), (1, 3, 2))

    #* merge state and output
    if isnothing(solved_states)
        final_output_arr = ele_output_arr
    else
        final_output_arr = cat(solved_states, ele_output_arr, dims=1)
    end
    final_output_arr
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
- `input`: Input data (Matrix for single node, Array for multiple nodes)
- `pas::ComponentVector`: Parameters and initial states
- `timeidx::Vector`: Time index vector
- `ptypes::Vector{Symbol}`: Parameter types (for multiple node input)
- `solver::AbstractSolver`: ODE solver to use (default: ODESolver())

# Returns
- `sol`: Solution of the ordinary differential equations

# Dimensions
- Input: var_names × node_names × ts_len (or var_names × ts_len for single node)
- Output: var_names × node_names × ts_len (or var_names × ts_len for single node)
"""
function solve_single_prob(
    ele::HydroBucket;
    input::Matrix,
    pas::ComponentVector,
    timeidx::Vector=collect(1:length(input[1])),
    solver::AbstractSolver=ODESolver(),
)
    params, init_states = pas[:params], pas[:initstates]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_list = map(eachrow(input)) do var
        LinearInterpolation(var, timeidx, extrapolate=true)
    end
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params)[1][nm].idx for nm in ele.infos[:param]]
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    ode_param_func = (p) -> p[:params][params_idx]

    if !isempty(ele.infos[:nn])
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in ele.infos[:nn]]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function single_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        du[:] = ele.ode_func(ode_input, u, ode_params, nn_params)
    end

    #* Solve the problem using the solver wrapper
    sol = solver(single_ele_ode_func!, pas, collect(init_states[ele.infos[:state]]), timeidx)
    if sol == false
        sol = zeros(length(ele.infos[:state]), length(timeidx))
    end
    sol
end

function solve_multi_prob(
    ele::HydroBucket;
    input::Array,
    pas::ComponentVector,
    timeidx::Vector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
    solver::AbstractSolver=ODESolver()
)
    #* 针对多个相同的state function采用并行化计算,这样能够避免神经网络反复多次计算减少梯度计算反馈
    #* 同时将多组state function放到同一个ode function中,这种并行计算或能提高预测性能,
    #* 这样每个步长的输入维度就是:节点个数*输入变量数
    #* 当前只针对unit相同的同步求解:
    params, init_states = pas[:params], pas[:initstates]

    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itpfunc_vecs = [LinearInterpolation.(eachslice(input[:, i, :], dims=1), Ref(timeidx), extrapolate=true) for i in 1:size(input)[2]]
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]

    #* 准备初始状态
    init_states_vec = collect([collect(init_states[ptype][ele.infos[:state]]) for ptype in ptypes])
    init_states_matrix = reduce(hcat, init_states_vec)

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params[ptypes[1]])[1][nm].idx for nm in ele.infos[:param]]
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    ode_param_func = (p) -> [p[:params][ptype][params_idx] for ptype in ptypes]

    #* 准备神经网络的参数
    if length(ele.infos[:nn]) > 0
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in ele.infos[:nn]]
        ode_nn_param_func = (p) -> [p[:nn][idx] for idx in nn_params_idx]
    else
        ode_nn_param_func = (_) -> nothing
    end

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function multi_ele_ode_func!(du, u, p, t)
        ode_input = ode_input_func(t)
        ode_params = ode_param_func(p)
        nn_params = ode_nn_param_func(p)
        tmp_output_vec = ele.ode_func.(ode_input, eachslice(u, dims=2), ode_params, Ref(nn_params))
        tmp_output = reduce(hcat, tmp_output_vec)
        du[:] = tmp_output
    end

    #* Solve the problem using the solver wrapper
    sol = solver(multi_ele_ode_func!, pas, init_states_matrix, timeidx)
    if sol == false
        sol = zeros(length(ele.infos[:state]), length(timeidx), length(ptypes))
    end
    sol
end

