"""
    HydroBucket(name::Symbol; funcs::Vector, dfuncs::Vector=StateFlux[])

Represents a hydrological bucket model component.

# Arguments
- `name::Symbol`: A symbol representing the name of the HydroBucket instance. If not provided, a name will be automatically generated from state variable names.
- `funcs::Vector`: A vector of flux functions that describe the hydrological processes.
- `dfuncs::Vector`: A vector of state derivative functions (default is an empty vector of StateFlux).

# Fields
- `funcs::Vector{<:AbstractFlux}`: Vector of flux functions describing hydrological processes.
- `dfuncs::Vector{<:AbstractStateFlux}`: Vector of state derivative functions for ODE calculations.
- `flux_func::Function`: Combined function for calculating all hydrological fluxes.
- `ode_func::Union{Nothing,Function}`: Function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed.
- `meta::HydroMeta`: Contains metadata about the bucket, including input, output, state, parameter, and neural network names.

# Description
HydroBucket is a structure that encapsulates the behavior of a hydrological bucket model. 
It combines multiple flux functions and state derivative functions to model water movement 
and storage within a hydrological unit.

The structure automatically extracts relevant information from the provided functions to 
populate the metadata, which includes names of:
- Inputs: Variables that drive the model
- Outputs: Variables produced by the model
- States: Internal model states that evolve over time
- Parameters: Model parameters that control behavior
- Neural Networks: Any neural network components (if applicable)

The `flux_func` and `ode_func` are constructed based on the provided `funcs` and `dfuncs`, 
enabling efficient calculation of fluxes and state changes over time.

This structure is particularly useful for building complex hydrological models by combining 
multiple HydroBucket instances to represent different components of a water system.

"""
struct HydroBucket{S,N} <: AbstractBucket
    """
    Vector of flux functions describing hydrological processes.
    """
    funcs::Vector{<:AbstractHydroFlux}
    """
    Vector of state derivative functions for ODE calculations.
    """
    dfuncs::Vector{<:AbstractStateFlux}
    """
    Generated function for calculating all hydrological fluxes.
    """
    flux_func::Function
    """
    Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed.
    """
    ode_func::Union{Nothing,Function}
    """
    Metadata about the bucket, including input, output, state, parameter, and neural network names.
    """
    meta::HydroMeta

    function HydroBucket(;
        funcs::Vector{<:AbstractHydroFlux},
        dfuncs::Vector{<:AbstractStateFlux}=StateFlux[],
        name::Union{Symbol,Nothing}=nothing,
        sort_funcs::Bool=false
    )
        funcs = sort_funcs ? sort_funcs(funcs) : funcs
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(funcs, dfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(vcat(funcs, dfuncs))
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(funcs)
        #* Setup the name information of the hydrobucket
        bucket_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), state_names)), :_bucket) : name
        meta = HydroMeta(bucket_name, input_names, output_names, param_names, state_names, nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(funcs, dfuncs, meta)
        return new{length(state_names) > 0 ? Tuple(state_names) : nothing, length(nn_names) > 0}(
            funcs,
            dfuncs,
            flux_func,
            ode_func,
            meta,
        )
    end
end

"""
    (ele::HydroBucket)(input::Matrix, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)
    (ele::HydroBucket)(input::Array, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)
    (ele::HydroBucket)(input::NamedTuple, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)
    (ele::HydroBucket)(input::Vector{<:NamedTuple}, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Run the HydroBucket model for given input and parameters.

# Arguments
- `input`: Input data. Can be a Matrix (dims: variables × time) or Array (dims: variables × nodes × time)
           or NamedTuple (variables => time series) or Vector{NamedTuple} (variables => time series, for multiple nodes)
- `pas`: ComponentVector containing parameters and initial states
- `config`: Configuration options including:
  - `solver`: AbstractHydroSolver to use for ODE solving (default: ODESolver())
  - `interp`: Interpolation method for input data (default: LinearInterpolation)
  - `timeidx`: Vector of time indices (default: 1:size(input, last_dim))
  - `ptypes`: (Array input only) Parameter types to use (default: all parameter types)
  - `stypes`: (Array input only) State types to use (default: all state types)
- `kwargs`: Additional keyword arguments
  - `convert_to_ntp`: Whether to convert output to NamedTuple (default: false)

# Returns
If convert_to_ntp=false (default):
- Matrix (dims: (states+outputs) × time) for Matrix input
- Array (dims: (states+outputs) × nodes × time) for Array input

If convert_to_ntp=true:
- NamedTuple of time series for Matrix input
- Vector of NamedTuples for Array input

# Details
- For Matrix input: Processes single node/location data
- For Array input: Processes multiple nodes/locations simultaneously
- If the bucket has an ODE function, solves states over time
- Calculates fluxes using the model's flux function
- Concatenates solved states (if any) with calculated fluxes for output
- Input dimensions must match number of input variables defined in model
- Required parameters and initial states must be present in pas
"""
function (ele::HydroBucket{S,N})(
    input::AbstractArray{T,2},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
) where {S,N,T}
    #* get kwargs
    solver = get(config, :solver, ManualSolver{true}())
    interp = get(config, :interp, DataInterpolations.LinearInterpolation)
    timeidx = get(config, :timeidx, collect(1:size(input, 2)))

    initstates_vec = pas[:initstates]
    model_params = pas[:params]
    nn_params = ifelse(N, pas[:nns], Vector{eltype(pas)}[])

    #* prepare input interpolation
    itpfunc_list = map((var) -> interp(var, timeidx), eachrow(input))
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]

    #* prepare parameter and nn parameter
    pas_vec = vcat(nn_params, model_params)
    nn_idx_bounds = 1:length(nn_params)
    params_idx_bound = length(nn_params)+1:length(pas_vec)

    #* define the ODE function
    function du_func(u,p,t)
        @views ps, nn_ps = p[params_idx_bound], p[nn_idx_bounds]
        ele.ode_func(ode_input_func(t), u, ps, nn_ps)
    end

    #* solve the problem by call the solver
    solved_states = solver(du_func, pas_vec, initstates_vec, timeidx)
    #* calculate output, slice input on time dim, then calculate each output
    tmp_flux_func = (i,s) -> ele.flux_func(i, s, model_params, nn_params)
    flux_output = tmp_flux_func.(eachslice(input, dims=2), eachslice(solved_states, dims=2))
    #* convert vector{vector} to matrix
    reduce(hcat, flux_output)
end

function (ele::HydroBucket{nothing,N})(
    input::AbstractArray{T,2},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
) where {N,T}
    model_params = pas[:params]
    nn_params = ifelse(N, pas[:nns], Vector{eltype(pas)}[])
    #* calculate output, slice input on time dim, then calculate each output
    tmp_flux_func = (i) -> ele.flux_func(i, nothing, model_params, nn_params)
    flux_output = tmp_flux_func.(eachslice(input, dims=2))
    #* convert vector{vector} to matrix
    reduce(hcat, flux_output)
end

function (ele::HydroBucket{S,N})(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
) where {S,N,T}
    input_dims, num_nodes, time_len = size(input)
    #* get kwargs
    solver = get(config, :solver, ManualSolver{true}())
    interp = get(config, :interp, LinearInterpolation)
    timeidx = get(config, :timeidx, Vector(1:size(input, 3)))

    #* prepare parameter and nn parameter
    params_len = length(get_param_names(ele))
    initstates_mat = pas[:initstates]

    params_mat = pas[:params]
    nn_params_vec = if N
        vec(pas[:nns])
    else 
        Vector{eltype(pas)}[]
    end
    vcat_pas = vcat(nn_params_vec, vec(params_mat))
    nn_idx_bounds = 1:length(nn_params_vec)
    params_idx_bound = length(nn_params_vec)+1:length(vcat_pas)

    #* prepare input function
    itpfunc_vecs = [interp.(eachslice(input_, dims=1), Ref(timeidx), extrapolate=true) for input_ in eachslice(input, dims=2)]
    ode_input_func = (t) -> [[itpfunc(t) for itpfunc in itpfunc_vec] for itpfunc_vec in itpfunc_vecs]

    #* define the ODE function
    function du_func(u,p,t)
        @views ps, nn_ps = reshape(p[params_idx_bound], params_len, num_nodes), p[nn_idx_bounds]
        ele.ode_func.(ode_input_func(t), eachslice(u, dims=2), eachslice(ps, dims=2), Ref(nn_ps))
    end

    #* Call the solve_prob method to solve the state of bucket at the specified timeidx
    solved_states = solver(du_func, vcat_pas, initstates_mat, timeidx)

    #* run other functions
    tmp_flux_func(i,s,p) = ele.flux_func(i, s, p, nn_params)
    ele_output_vec = map(1:size(input, 3)) do i
        input_ = @view input[:, :, i]
        states_ = @view solved_states[:, :, i]
        reduce(hcat, ele.flux_func.(eachslice(input_, dims=2), eachslice(states_, dims=2), eachslice(params_mat, dims=2), Ref(nn_params_vec)))
    end
    reduce((m1, m2) -> cat(m1, m2, dims=3), ele_output_vec)
end

function (ele::HydroBucket{nothing,N})(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
) where {N,T}
    initstates_mat = pas[:initstates]
    params_mat = pas[:params]
    nn_params_vec = if N
        vec(pas[:nns])
    else 
        Vector{eltype(pas)}[]
    end
    #* run other functions
    tmp_flux_func(i,p) = ele.flux_func(i, nothing, p, nn_params)
    ele_output_vec = map(1:size(input, 3)) do i
        input_ = @view input[:, :, i]
        states_ = @view solved_states[:, :, i]
        reduce(hcat, ele.flux_func.(eachslice(input_, dims=2), eachslice(states_, dims=2), eachslice(params_mat, dims=2), Ref(nn_params_vec)))
    end
    reduce((m1, m2) -> cat(m1, m2, dims=3), ele_output_vec)
end