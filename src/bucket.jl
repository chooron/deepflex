"""
$(TYPEDEF)
The basic hydrological calculation module contains multiple hydrological fluxes,
and can simulate the balance calculation of a physical module.
# Fields
$(FIELDS)
# Example
```
```
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
    input::AbstractMatrix,
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
    flux_output_matrix = run_fluxes(ele, input=fluxes, pas=pas)
    #* merge output and state
    if isnothing(solved_states)
        output_matrix = flux_output_matrix
    else
        output_matrix = cat(solved_states, flux_output_matrix, dims=1)
    end
    output_matrix
end

function run_fluxes(
    ele::HydroBucket;
    #* var num * ts len
    input::AbstractMatrix,
    pas::ComponentVector,
)
    params_vec = collect([pas[:params][nm] for nm in ele.infos[:param]])
    if length(ele.infos[:nn]) > 0
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.infos[:nn]])
    else
        nn_params_vec = nothing
    end
    ele_output = ele.flux_func.(eachcol(input), Ref(params_vec), Ref(nn_params_vec))
    #* convert vector{vector} to matrix
    ele_output_matrix = reduce(hcat, ele_output)
    ele_output_matrix
end

function solve_single_prob(
    ele::HydroBucket;
    input::AbstractMatrix,
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
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_param_func = (p) -> [p[:params][idx] for idx in params_idx]

    if length(ele.infos[:nn]) > 0
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
    sol
end

"""
input dims: var_names * node_names * ts_len
pas type: ComponentVector(params=(node_1=(p1=,p2=,)), initstates=(...), nn=(...))

output dims: var_names * node_names * ts_len
"""
function run_multi_fluxes(
    ele::HydroBucket;
    #* var num * node num * ts len
    input::AbstractArray,
    pas::ComponentVector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
)
    ele_param_vec = collect([collect([pas[:params][ptype][pname] for pname in ele.infos[:param]]) for ptype in ptypes])
    if length(ele.infos[:nn]) > 0
        nn_params_vec = collect([pas[:nn][nm] for nm in ele.infos[:nn]])
    else
        nn_params_vec = nothing
    end
    #* array dims: (num of node, sequence length, variable dim)
    ele_output_vec = [ele.flux_func.(eachslice(input[:, i, :], dims=2), Ref(ele_param_vec[i]), Ref(nn_params_vec)) for i in 1:size(input)[2]]
    ele_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, u) for u in ele_output_vec])
    permutedims(ele_output_arr, (1, 3, 2))
end

"""
input dims: var_names * node_names * ts_len
pas type: ComponentVector(params=(node_1=(p1=,p2=,)), initstates=(...), nn=(...))

output dims: var_names * node_names * ts_len
"""
function solve_multi_prob(
    ele::HydroBucket;
    input::AbstractArray,
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
    itpfunc_vecs = [LinearInterpolation.(eachslice(i, dims=1), Ref(timeidx), extrapolate=true) for i in eachslice(input, dims=2)]
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]

    #* 准备初始状态
    init_states_vec = collect([collect(init_states[ptype][ele.infos[:state]]) for ptype in ptypes])
    init_states_matrix = reduce(hcat, init_states_vec)

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    params_idx = [getaxes(params[ptypes[1]])[1][nm].idx for nm in ele.infos[:param]]
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func (see in line 45)
    ode_param_func = (p) -> [[p[:params][ptype][idx] for idx in params_idx] for ptype in ptypes]

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
    sol
end

