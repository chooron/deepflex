"""
    HydroFlux

Represents a simple flux component in a hydrological model.

# Arguments
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and parameter names.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and parameter names.

# Constructors
    HydroFlux(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num}, meta::HydroMeta)
    HydroFlux(flux_names::Pair{Vector{Symbol},Vector{Symbol}}, param_names::Vector{Symbol}=Symbol[]; flux_funcs::Vector{<:Function}=Function[])

# Description
`HydroFlux` is a structure that encapsulates a simple flux calculation in a hydrological model. 
It can be constructed either by providing explicit inputs, outputs, parameters, and expressions, 
or by specifying names for fluxes and parameters along with optional flux functions.

The structure automatically compiles the provided expressions or functions into an efficient 
calculation function, which can be used to compute flux values given input and parameter values.

This structure is particularly useful for representing straightforward hydrological processes 
where the relationship between inputs and outputs can be expressed as simple mathematical formulas.
"""
struct HydroFlux{T<:Num,F<:Function,M<:HydroMeta} <: AbstractHydroFlux
    "Vector of input variables"
    inputs::Vector{T}
    "Vector of output variables"
    outputs::Vector{T}
    "Vector of parameter variables"
    params::Vector{T}
    "Vector of expressions describing the formulas for output variables"
    exprs::Vector{T}
    "Compiled function that calculates the flux"
    func::F
    "Metadata about the flux, including input, output, and parameter names"
    meta::M

    function HydroFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        params::Vector{T};
        exprs::Vector{T}=T[],
    ) where {T<:Num}
        #* name the flux
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        flux_name = Symbol(Symbol(reduce((x, y) -> Symbol(x, y), output_names)), :_sflux)
        #* construct meta
        meta = HydroMeta(name=flux_name, inputs=inputs, outputs=outputs, params=params)
        #* if no expression provided, use the hydrology formula library to build the flux
        if length(exprs) == 0
            @info "No expression provided, using the hydrology formula library (`HydroModelLibrary.jl`) to build the flux"
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydrofunc = HydroEquation(input_names, output_names, param_names)
            exprs = HydroModelLibrary.expr(hydrofunc)
        end
        #* build flux function
        flux_func = build_flux_func(inputs, outputs, params, exprs)

        return new{T,typeof(flux_func),typeof(meta)}(inputs, outputs, params, exprs, flux_func, meta)
    end

    #* construct hydro flux with input fluxes and output fluxes
    HydroFlux(fluxes::Pair{Vector{Num},Vector{Num}}, params::Vector{Num}=Num[]; exprs::Vector{Num}=Num[]) = HydroFlux(fluxes[1], fluxes[2], params, exprs=exprs)
end

"""
    (flux::AbstractHydroFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the simple flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API
    - `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).
"""
function (flux::AbstractHydroFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, 1)
    params_vec = collect([pas[:params][nm] for nm in get_param_names(flux)])
    flux.func(input, params_vec, timeidx)
end

function (flux::AbstractHydroFlux)(input::AbstractArray{N,2}, pas::ComponentVector; kwargs...) where {N}
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    # assert the input params must include all the parameters in the flux
    @assert length(timeidx) == size(input, 2) "Time index length does not match the number of time steps"
    @assert all(nm in keys(pas[:params]) for nm in get_param_names(flux)) "Input parameters do not match the flux parameters, the flux parameters should be: $(get_param_names(flux))"
    params_vec = collect([pas[:params][nm] for nm in get_param_names(flux)])
    output_arr = reduce(hcat, flux.func.(eachslice(input, dims=2), Ref(params_vec), timeidx))
    output_arr
end

function (flux::AbstractHydroFlux)(input::AbstractArray{N,3}, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...) where {N}
    #* get kwargs
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))
    @assert length(timeidx) == size(input, 3) "Time index length does not match the number of time steps"

    #* extract params and nn params
    params_vec = collect([collect([pas[:params][ptype][pname] for pname in get_param_names(flux)]) for ptype in ptypes])

    #* array dims: (var_names * node_names * ts_len)
    flux_output_vec = [reduce(hcat, flux.func.(eachslice(input[:, :, i], dims=2), params_vec, timeidx[i])) for i in eachindex(timeidx)]
    if length(flux_output_vec) == 1
        tmp_output_arr = flux_output_vec[1]
        flux_output_arr = reshape(tmp_output_arr, size(tmp_output_arr)..., 1)
    else
        flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    end
    flux_output_arr
end


"""
    NeuralFlux <: AbstractNeuralFlux

Represents a neural network-based flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `nnparam::Symbolics.Arr`: An array of neural network parameters.
- `expr::Symbolics.Arr{Num,1}`: Expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux using the neural network.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and neural network parameter names.
- `nninfos::NamedTuple`: Contains information about the neural network's input and output structure.

# Constructors
    # 1. Construct a NeuralFlux with specified input/output fluxes and a neural network, the neural network should specify the name
    NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain::Lux.AbstractExplicitContainerLayer)

# Description
`NeuralFlux` is a structure that encapsulates a neural network-based flux calculation in a hydrological model. 
It combines symbolic mathematics with neural networks to create a flexible and powerful representation of complex hydrological processes.

The structure automatically handles the integration of the neural network into the symbolic framework, 
allowing for seamless incorporation of machine learning models into traditional hydrological equations.

This structure is particularly useful for representing complex, non-linear relationships in hydrological systems 
where traditional equations may be insufficient or unknown.
"""
struct NeuralFlux{T<:Num,F<:Function,M<:HydroMeta} <: AbstractNeuralFlux
    "Vector of input variables"
    inputs::Vector{T}
    "Vector of output variables"
    outputs::Vector{T}
    "Array of neural network parameters"
    nnparam::Symbolics.Arr{T,1}
    "Array of expressions describing the formulas for output variables"
    expr::Symbolics.Arr{T,1}
    "Compiled function that calculates the flux using the neural network"
    func::F
    "Metadata about the flux, including input, output, and neural network parameter names"
    meta::M
    "Information about the neural network's input and output structure"
    nninfos::NamedTuple

    function NeuralFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        chain::LuxCore.AbstractLuxLayer
    ) where {T<:Num}
        #* assert the chain has a name
        @assert chain.name isa Symbol "Neural network chain should have a name with Symbol type"
        #* Get the neural network name (neural flux param name) and object
        chain_name = chain.name
        #* Initialize model parameter type for model parameter dimension definition
        init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
        params_axes = getaxes(init_params)

        #* Define parameter variables according to initialization parameters: Define type as Vector{parameter length}
        chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(init_params))
        #* Use Symbolics.array_term to define the slow-building parameter variables:
        #* when the model is called, it is rebuilt into the ComponentVector type according to
        #* the axes of `init_params` and the Vector type of the parameter as the calculation parameter input
        lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, params_axes, size=size(chain_params))

        #* Constructing neural network input and output variables
        #* The input and output of the model can only be of type Symbolics.Arr{Num, 1},
        #* so it cannot be defined based on outputs and output_vars
        nn_input_name = Symbol(chain_name, :_input)
        nn_output_name = Symbol(chain_name, :_output)
        nn_input = first(@variables $(nn_input_name)[1:length(inputs)])
        nn_output = first(@variables $(nn_output_name)[1:length(outputs)])

        #* Constructing a calculation expression based on a neural network
        flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)
        nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, params_axes))

        #* neuralflux meta
        meta = HydroMeta(name=Symbol(chain_name, :_nflux), inputs=inputs, outputs=outputs, nn_names=[chain_name])
        nninfos = (inputs=nn_input, outputs=nn_output, paramlen=length(init_params))
        new{T,typeof(nn_func),typeof(meta)}(
            inputs, outputs, chain_params,
            flux_expr, nn_func,
            meta, nninfos,
        )
    end

    function NeuralFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        params::Symbolics.Arr{T,1},
        expr::Symbolics.Arr{T,1},
        chain,
    ) where {T<:Number}
        input_vec, output_vec = Symbol(chain.name, :_input), Symbol(chain.name, :_output)
        nninfos = (inputs=first(@variables $(input_vec)[1:length(inputs)]), outputs=first(@variables $(output_vec)[1:length(outputs)]))
        meta = HydroMeta(name=Symbol(chain.name, :_nflux), inputs=inputs, outputs=outputs, nn_names=[chain.name])
        params_axes = getaxes(ComponentVector(Lux.initialparameters(StableRNG(42), chain)))
        nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, params_axes))
        new{T,typeof(nn_func),typeof(meta)}(inputs, outputs, params, expr, nn_func, meta, nninfos)
    end

    #* construct neural flux with input fluxes and output fluxes
    NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain) = NeuralFlux(fluxes[1], fluxes[2], chain)
end

"""
    (flux::AbstractFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the flux model (simple or neural) to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).

# Note
For neural flux models, the parameters are accessed from `pas[:nn]` instead of `pas[:params]`.
"""
function (flux::AbstractNeuralFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][get_nn_names(flux)[1]]
    flux.func(input, nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T}
    nn_params_vec = pas[:nn][get_nn_names(flux)[1]]
    output_arr = flux.func(input, nn_params_vec)
    output_arr
end

function (flux::AbstractNeuralFlux)(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...) where {T}
    nn_params = pas[:nn][get_nn_names(flux)[1]]
    #* array dims: (ts_len * node_names * var_names)
    flux_output_vec = [flux.func(input[:, i, :], nn_params) for i in 1:size(input)[2]]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (1, 3, 2))
end

"""
    StateFlux <: AbstractStateFlux

Represents a state flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `state::Num`: The state variable.
- `params::Vector{Num}`: A vector of parameter variables.
- `expr::Num`: The expression describing the state variable's formula.
- `func::Function`: A function to calculate the state flux.
- `meta::HydroMeta`: Contains metadata about inputs, state, parameters, and neural networks.

# Constructors
    # 1. Detailed specification of inputs, state, parameters, and state expression
    StateFlux(fluxes::Vector{Num}, state::Num, params::Vector{Num}=Num[]; expr::Num)
    # 2. Automatic construction of state expression as the difference between sum of input fluxes and sum of output fluxes
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num)

# Description
StateFlux is a structure that represents a state flux in a hydrological model. It encapsulates 
the relationship between input fluxes, output fluxes, and a state variable. The structure 
provides methods to calculate state changes based on the provided expressions and parameters.

The first constructor allows for detailed specification of inputs, state, parameters, and the 
state expression. The second constructor automatically constructs a state expression as the 
difference between sum of input fluxes and sum of output fluxes. The third constructor is used 
for simple state transitions.

This structure is particularly useful in building complex hydrological models where state 
variables evolve over time based on various input and output fluxes.
"""
struct StateFlux{T<:Num,F<:Function,M<:HydroMeta} <: AbstractStateFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{T}
    "A map of state names (Symbol) and its variables (Num)"
    state::T
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{T}
    "flux expressions to descripe the formula of the state variable"
    expr::T
    "flux expressions to descripe the formula of the output variable"
    func::F
    "bucket information: keys contains: input, output, param, state"
    meta::M

    function StateFlux(
        inputs::Vector{T}, state::T, params::Vector{T}=T[]; expr::T
    ) where {T<:Num}
        #* Convert to a symbol based on the variable
        state_name = Symbolics.tosymbol(state, escape=false)
        meta = HydroMeta(name=Symbol(state_name, :_stflux), inputs=inputs, states=[state], params=params)
        state_func = build_flux_func(inputs, [state], params, [expr])
        return new{T,typeof(state_func),typeof(meta)}(
            inputs,
            state,
            params,
            expr,
            state_func,
            meta
        )
    end
    #* construct state flux with input fluxes and output fluxes
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num) = StateFlux(vcat(fluxes[1], fluxes[2]), state, expr=sum(fluxes[1]) - sum(fluxes[2]))
    #* construct state flux with state variables
    StateFlux(states::Pair{Num,Num}) = StateFlux([states[2]], states[1], expr=states[2] - states[1])
end

"""
    (flux::AbstractStateFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the state flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
This function does not actually return a value, as state flux models cannot be run directly.

# Notes
- State flux models cannot be run directly and will throw an error if attempted.
- To use state flux models, they should be incorporated into a HydroFlux or other composite flux model.
"""
(::AbstractStateFlux)(::AbstractVector, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::AbstractArray{T,2}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::AbstractArray{T,3}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"

