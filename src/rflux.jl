
"""
    RouteFlux{rtype} <: AbstractRouteFlux

Represents a routing flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `meta::NamedTuple`: Contains metadata about the flux, including input, output, and parameter names.

# Type Parameters
- `rtype`: A symbol representing the specific type of routing flux.

# Constructors
    RouteFlux(input::Num, params::Vector{Num}; routetype::Symbol, output::Union{Num,Nothing}=nothing)

# Description
`RouteFlux` is a structure that encapsulates a routing flux calculation in a hydrological model. 
It is designed to represent various types of routing processes, with the specific type indicated 
by the `rtype` parameter.

The structure automatically generates an output variable name based on the input if not provided,
and organizes the information about inputs, outputs, and parameters into the `meta` field.

This structure is particularly useful for representing routing processes in hydrological models, 
where water is transferred from one point to another in the system.
"""
struct StateRouteFlux{rtype} <: AbstractStateRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "A map of states names (Symbol) and its variables (Num)"
    states::Vector{Num}
    "bucket information: keys contains: input, output, param, state"
    meta::HydroMeta

    function StateRouteFlux(
        input::Num,
        params::Vector{Num},
        states::Vector{Num};
        routetype::Symbol,
        output::Union{Num,Nothing}=nothing,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_names = isempty(params) ? Symbol[] : Symbolics.tosymbol.(params, escape=false)
        state_names = isempty(states) ? Symbol[] : Symbolics.tosymbol.(states, escape=false)

        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end
        #* Setup the name information of the hydroroutement
        meta = HydroMeta(inputs=[input_name], outputs=[output_name], params=param_names, states=state_names, name=Symbol(output_name, :_route_flux))

        return new{routetype}(
            [input],
            [output],
            params,
            states,
            meta
        )
    end
end

function DischargeRouteFlux(input::Num, output::Union{Num,Nothing}=nothing)
    @parameters lag
    @variables s_river

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return StateRouteFlux(input, [lag], [s_river]; routetype=:discharge, output=output)
end

function RiverRouteFlux(input::Num, output::Union{Num,Nothing}=nothing)
    @parameters lag
    @variables s_river

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return StateRouteFlux(input, [lag], [s_river]; routetype=:river, output=output)
end

struct DynamicRouteFlux{rtype} <: AbstractDynamicRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "bucket information: keys contains: input, output, param, state"
    meta::HydroMeta

    function DynamicRouteFlux(
        input::Num,
        params::Vector{Num};
        routetype::Symbol,
        output::Union{Num,Nothing}=nothing,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_names = isempty(params) ? Symbol[] : Symbolics.tosymbol.(params, escape=false)

        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end
        #* Setup the name information of the hydroroutement
        meta = HydroMeta(inputs=[input_name], outputs=[output_name], params=param_names, name=Symbol(output_name, :_route_flux))

        return new{routetype}(
            [input], [output],
            params,
            meta
        )
    end
end

function MuskingumRouteFlux(input::Num, output::Union{Num,Nothing}=nothing)
    @parameters k x

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return DynamicRouteFlux(input, [k, x]; routetype=:muskingum, output=output)
end


"""
    get_rflux_func(::AbstractRouteFlux; pas::ComponentVector, ptypes::AbstractVector{Symbol})

Get the function for calculating the outflow process of a route flux model.

This function should be implemented by different types of RouteFlux. It returns a function that calculates the outflow process.

The returned function should have the following input parameters:
- `du`: Change in state variables
- `s_rivers`: State values of the RouteFlux
- `q_in`: Input flow for each time step
- `q_gen`: Generated runoff for each time step (after area conversion)
- `p`: Parameters

The returned function should output:
- `q_out`: Outflow for each time step
- Updates to `du` representing changes in `s_rivers`

# Arguments
- `::AbstractRouteFlux`: The RouteFlux instance
- `pas::ComponentVector`: Component vector containing parameters
- `ptypes::AbstractVector{Symbol}`: Vector of parameter type symbols

# Returns
A function that calculates the outflow process based on the specific RouteFlux implementation.

# Note
This function must be implemented by subtypes of AbstractRouteFlux to provide the specific outflow calculation logic for each routing method.
"""
get_rflux_func(::AbstractRouteFlux; kwargs...) = @error "Must be implemented by subtype"

function get_rflux_func(::StateRouteFlux{:discharge}; kwargs...)

    function rflux_func(s_river, q_in, q_gen, p)
        q_rf = @. (s_river + q_in) / (1 + p.lag)
        d_state = q_in .- q_rf
        q_out = q_rf .+ q_gen
        q_out, d_state
    end

    return rflux_func
end

function get_rflux_func(::StateRouteFlux{:river}; kwargs...)

    function rflux_func(s_river, q_in, q_gen, p)
        k_ps, x_ps = p.k, p.x
        q_rf = @.((s_river - k_ps * x_ps * q_in) / (k_ps * (1 - x_ps)))
        d_state = q_in .- q_rf
        q_out = q_rf .+ q_gen
        q_out, d_state
    end

    return rflux_func
end

function get_rflux_func(::DynamicRouteFlux{:muskingum}; adjacency::AbstractMatrix, kwargs...)

    function rflux_func(q_out_t1, q_gen, p)
        c0, c1, c2 = p.c0, p.c1, p.c2
        b = c0 * q_gen + c1 * (adjacency * q_out_t1 + q_gen) + c2 * q_out_t1
        return b
    end

    return rflux_func
end

"""
    get_rflux_initstates(::AbstractRouteFlux; input::AbstractMatrix, pas::ComponentVector, ptypes::AbstractVector{Symbol}, stypes::AbstractVector{Symbol})

Get the initial states for a route flux model.

This function is designed to accommodate different subclasses of `AbstractRouteFlux` that may have varying methods of constructing initial states. It provides a flexible interface for generating initial state vectors based on input data, parameters, parameter types, and state types.

# Arguments
- `::AbstractRouteFlux`: The route flux model instance.
- `input::AbstractMatrix`: Input data matrix, typically representing time series of hydrological variables.
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter types.
- `stypes::AbstractVector{Symbol}`: A vector of symbols representing state types.

# Returns
A vector representing the initial states for the route flux model.

# Notes
- This function must be implemented by subtypes of `AbstractRouteFlux` to provide specific initialization logic for each routing method.
- The function allows for flexibility in how initial states are constructed, which can vary depending on the specific requirements of each routing method.
- By taking input data, parameters, and type information as arguments, the function can create initial states that are appropriately tailored to the specific model configuration and data characteristics.

"""
get_rflux_initstates(::AbstractStateRouteFlux; kwargs...) = @error "Must be implemented by subtype"
function get_rflux_initstates(::StateRouteFlux{:discharge}; pas::ComponentVector, stypes::AbstractVector{Symbol}, kwargs...)
    [pas[:initstates][stype][:s_river] for stype in stypes]
end
function get_rflux_initstates(::StateRouteFlux{:river}; pas::ComponentVector, stypes::AbstractVector{Symbol}, kwargs...)
    [pas[:initstates][stype][:s_river] for stype in stypes]
end
"""
    get_rflux_initstates(::AbstractRouteFlux; input::AbstractMatrix, pas::ComponentVector, ptypes::AbstractVector{Symbol}, stypes::AbstractVector{Symbol})

Get the initial states for a route flux model.

This function is designed to accommodate different subclasses of `AbstractRouteFlux` that may have varying methods of constructing initial states. It provides a flexible interface for generating initial state vectors based on input data, parameters, parameter types, and state types.

# Arguments
- `::AbstractRouteFlux`: The route flux model instance.
- `input::AbstractMatrix`: Input data matrix, typically representing time series of hydrological variables.
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter types.
- `stypes::AbstractVector{Symbol}`: A vector of symbols representing state types.

# Returns
A vector representing the initial states for the route flux model.

# Notes
- This function must be implemented by subtypes of `AbstractRouteFlux` to provide specific initialization logic for each routing method.
- The function allows for flexibility in how initial states are constructed, which can vary depending on the specific requirements of each routing method.
- By taking input data, parameters, and type information as arguments, the function can create initial states that are appropriately tailored to the specific model configuration and data characteristics.

"""
get_rflux_parameters(::AbstractRouteFlux; kwargs...) = @error "Must be implemented by subtype"
function get_rflux_parameters(::StateRouteFlux{:discharge}; pas::ComponentVector, ptypes::AbstractVector{Symbol}, kwargs...)
    lag_ps = [pas[:params][ptype][:lag] for ptype in ptypes]
    return ComponentVector(lag=lag_ps)
end
function get_rflux_parameters(::StateRouteFlux{:river}; pas::ComponentVector, ptypes::AbstractVector{Symbol}, kwargs...)
    k_ps = [pas[:params][ptype][:k] for ptype in ptypes]
    x_ps = [pas[:params][ptype][:x] for ptype in ptypes]
    return ComponentVector(k=k_ps, x=x_ps)
end
function get_rflux_parameters(::DynamicRouteFlux{:muskingum}; pas::ComponentVector, ptypes::AbstractVector{Symbol}, delta_t::Number, adjacency::AbstractMatrix, kwargs...)
    k_ps = [pas[:params][ptype][:k] for ptype in ptypes]
    x_ps = [pas[:params][ptype][:x] for ptype in ptypes]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    A = Matrix(I, size(adjacency)...) .- diagm(c0) * adjacency
    return ComponentVector(c0=diagm(c0), c1=diagm(c1), c2=diagm(c2), A=A)
end

"""
    (flux::AbstractRouteFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the route flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
This function does not actually return a value, as route flux models are abstract and cannot be run directly.

# Notes
- Route flux models are abstract and cannot be run directly. Attempting to do so will throw an error.
- Specific implementations of route flux models (subtypes of AbstractRouteFlux) should provide their own implementations of this function.
- Route flux models are typically used to represent water routing processes in hydrological systems.
"""
(::AbstractRouteFlux)(::Vector, ::ComponentVector; kwargs...) = @error "Abstract RouteFlux is not support for single timepoint"
(::AbstractRouteFlux)(input::Matrix, pas::ComponentVector; kwargs...) = @error "Must be implemented by subtype"
(::AbstractRouteFlux)(input::Array, pas::ComponentVector; kwargs...) = @error "Must be implemented with the Route type"

# todo 直接调用还不完善
function (flux::StateRouteFlux{:discharge})(input::Matrix, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, collect(1:size(input)[2]))
    delta_t = get(kwargs, :delta_t, 1.0)
    solver = get(kwargs, :solver, ODESolver(alg=Rosenbrock23()))
    converted_timeidx = timeidx .* delta_t
    input_itp = LinearInterpolation(input[1, :], converted_timeidx)
    params = pas[:params]

    function msk_func!(du, u, p, t)
        s_river = u[1]
        q_in = input_itp(t)
        k, x = p[:params]
        q_out = (s_river - k * x * q_in) / (k * (1 - x))
        du[1] = q_in - q_out
    end

    init_states = [params.k * input[1, 1]]
    sol_arr = solver(msk_func!, pas, init_states, converted_timeidx)
    s_river_vec = sol_arr[1, :]
    q_out_vec = @.((s_river_vec - params.k * params.x * input[1, :]) / (params.k * (1 - params.x)))
    q_out_vec
end

function (flux::StateRouteFlux{:river})(input::Matrix, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, collect(1:size(input)[2]))
    delta_t = get(kwargs, :delta_t, 1.0)
    solver = get(kwargs, :solver, ODESolver(alg=Rosenbrock23()))
    converted_timeidx = timeidx .* delta_t
    input_itp = LinearInterpolation(input[1, :], converted_timeidx)
    params = pas[:params]

    function msk_func!(du, u, p, t)
        s_river = u[1]
        q_in = input_itp(t)
        k, x = p[:params]
        q_out = (s_river - k * x * q_in) / (k * (1 - x))
        du[1] = q_in - q_out
    end

    init_states = [params.k * input[1, 1]]
    sol_arr = solver(msk_func!, pas, init_states, converted_timeidx)
    s_river_vec = sol_arr[1, :]
    q_out_vec = @.((s_river_vec - params.k * params.x * input[1, :]) / (params.k * (1 - params.x)))
    q_out_vec
end

function (flux::DynamicRouteFlux{:muskingum})(input::Matrix, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, collect(1:size(input)[2]))
    delta_t = get(kwargs, :delta_t, 1.0)
    solver = get(kwargs, :solver, DiscreteSolver(alg=FunctionMap{true}()))
    input_itp = LinearInterpolation(input[1, :], timeidx .* delta_t)
    k, x = pas[:params]
    c0 = ((delta_t / k) - (2 * x)) / ((2 * (1 - x)) + (delta_t / k))
    c1 = ((delta_t / k) + (2 * x)) / ((2 * (1 - x)) + (delta_t / k))
    c2 = ((2 * (1 - x)) - (delta_t / k)) / ((2 * (1 - x)) + (delta_t / k))
    msk_params = ComponentVector(c0=c0, c1=c1, c2=c2)

    function msk_func!(du, u, p, t)
        c0, c1, c2 = p
        q_in_t1, q_out_t1 = u[1], u[2]
        q_in_t2 = input_itp(t)
        q_out_t2 = @. (c0 * q_in_t2) + (c1 * q_in_t1) + (c2 * q_out_t1)
        du[1] = q_in_t2 - q_in_t1
        du[2] = q_out_t2 - q_out_t1
    end

    init_states = [input[1, 1], input[1, 1]]
    sol_arr = solver(msk_func!, msk_params, init_states, timeidx)
    q_out_t2_vec = sol_arr[2, :]
    q_out_t2_vec
end