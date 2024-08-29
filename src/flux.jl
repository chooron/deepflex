"""
    SimpleFlux

Represents a simple flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux.
- `infos::NamedTuple`: Contains metadata about the flux, including input, output, and parameter names.

# Constructors
    SimpleFlux(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num}, infos::NamedTuple)
    SimpleFlux(flux_names::Pair{Vector{Symbol},Vector{Symbol}}, param_names::Vector{Symbol}=Symbol[]; flux_funcs::Vector{<:Function}=Function[])

# Description
`SimpleFlux` is a structure that encapsulates a simple flux calculation in a hydrological model. 
It can be constructed either by providing explicit inputs, outputs, parameters, and expressions, 
or by specifying names for fluxes and parameters along with optional flux functions.

The structure automatically compiles the provided expressions or functions into an efficient 
calculation function, which can be used to compute flux values given input and parameter values.

This structure is particularly useful for representing straightforward hydrological processes 
where the relationship between inputs and outputs can be expressed as simple mathematical formulas.
"""
struct SimpleFlux <: AbstractSimpleFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "flux expressions to descripe the formula of the output variable"
    exprs::Vector{Num}
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function SimpleFlux(
        inputs::Vector{Num},
        outputs::Vector{Num},
        params::Vector{Num},
        exprs::Vector{Num},
        infos::NamedTuple
    )
        flux_func = build_flux_func(inputs, outputs, params, exprs)

        return new(
            inputs,
            outputs,
            params,
            exprs,
            flux_func,
            infos
        )
    end

    function SimpleFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        param_names::Vector{Symbol}=Symbol[];
        flux_funcs::Vector{<:Function}=Function[],
    )
        #* Get input and output names
        input_names, output_names = flux_names[1], flux_names[2]
        infos = (input=input_names, output=output_names, param=param_names, nn=Symbol[], state=Symbol[])
        if length(flux_funcs) > 0
            #* Create variables by names
            inputs = [first(@variables $var) for var in input_names]
            outputs = [first(@variables $var) for var in output_names]
            params = [first(@parameters $var) for var in param_names]
            #* When a calculation function is provided, exprs are constructed based on the calculation function and variables
            exprs = [flux_func(inputs, params) for flux_func in flux_funcs]
        else
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            inputs, outputs, params = hydro_equation.inputs, hydro_equation.outputs, hydro_equation.params
            exprs = HydroEquations.expr(hydro_equation)
        end

        #* Building the struct
        return SimpleFlux(
            inputs,
            outputs,
            params,
            exprs,
            infos
        )
    end

    function SimpleFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        exprs::Vector{Num}=Num[]
    )
        #* Get input and output variables
        inputs, outputs = fluxes[1], fluxes[2]

        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = length(params) > 0 ? Symbolics.tosymbol.(params) : Symbol[]
        infos = (input=input_names, output=output_names, param=param_names, nn=Symbol[], state=Symbol[])

        if length(exprs) == 0
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            exprs = HydroEquations.expr(hydro_equation)
        end

        return SimpleFlux(
            inputs,
            outputs,
            params,
            exprs,
            infos
        )
    end
end

"""
    (flux::AbstractSimpleFlux)(input::Vector, pas::ComponentVector; kwargs...)

Apply the simple flux model to a vector input.

# Arguments
- `input::Vector`: A vector of input values.
- `pas::ComponentVector`: A component vector containing parameter values.
- `kwargs...`: Additional keyword arguments (unused).

# Returns
The result of applying the flux function to the input and parameters.
"""

"""
    (flux::AbstractSimpleFlux)(input::Matrix, pas::ComponentVector; kwargs...)

Apply the simple flux model to a matrix input.

# Arguments
- `input::Matrix`: A matrix of input values, where each column represents a different input vector.
- `pas::ComponentVector`: A component vector containing parameter values.
- `kwargs...`: Additional keyword arguments (unused).

# Returns
A matrix where each column is the result of applying the flux function to the corresponding input column.
"""

"""
    (flux::AbstractSimpleFlux)(input::Array, pas::ComponentVector, ptypes::AbstractVector{Symbol})

Apply the simple flux model to a 3D array input with multiple parameter types.

# Arguments
- `input::Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories.
- `kwargs...`: Additional keyword arguments (unused).

# Returns
A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).
"""
function (flux::AbstractSimpleFlux)(input::Vector, pas::ComponentVector; kwargs...)
    flux.func(input, collect([pas[:params][nm] for nm in flux.infos[:param]]), nothing)
end

function (flux::AbstractSimpleFlux)(input::Matrix, pas::ComponentVector; kwargs...)
    params_vec = collect([pas[:params][nm] for nm in flux.infos[:param]])
    reduce(hcat, flux.func.(eachslice(input, dims=2), Ref(params_vec), Ref(nothing)))
end

function (flux::AbstractSimpleFlux)(input::Array, pas::ComponentVector, ptypes::AbstractVector{Symbol}; kwargs...)
    param_vec = collect([collect([pas[:params][ptype][pname] for pname in flux.infos[:param]]) for ptype in ptypes])
    #* array dims: (var_names * node_names * ts_len)
    flux_output_vec = [reduce(hcat, flux.func.(eachslice(input[:, i, :], dims=2), Ref(param_vec[i]), Ref(nothing))) for i in eachindex(ptypes)]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (3, 1, 2))
end

"""
    StateFlux

Represents a state flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `state::Num`: The state variable.
- `params::Vector{Num}`: A vector of parameter variables.
- `expr::Num`: The expression describing the state variable's formula.
- `func::Function`: A function to calculate the state flux.
- `infos::NamedTuple`: Contains metadata about inputs, state, parameters, and neural networks.

# Constructors
    StateFlux(fluxes::Vector{Num}, state::Num, params::Vector{Num}=Num[]; expr::Num)
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num)
    StateFlux(states::Pair{Num,Num})

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

struct StateFlux <: AbstractStateFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of state names (Symbol) and its variables (Num)"
    state::Num
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "flux expressions to descripe the formula of the state variable"
    expr::Num
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function StateFlux(
        fluxes::Vector{Num},
        state::Num,
        params::Vector{Num}=Num[];
        expr::Num
    )
        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(fluxes, escape=false)
        state_name = Symbolics.tosymbol(state, escape=false)
        param_names = length(params) > 0 ? Symbol.(Symbolics.tosymbol.(params, escape=false)) : Symbol[]
        infos = (input=input_names, state=[state_name], param=param_names, nn=Symbol[])
        state_func = build_flux_func(fluxes, [state], params, [expr])
        return new(
            fluxes,
            state,
            params,
            expr,
            state_func,
            infos
        )
    end

    function StateFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        state::Num;
    )
        influxes, outfluxes = fluxes[1], fluxes[2]
        #* Construct the default calculation formula: sum of input variables minus sum of output variables
        state_expr = sum(influxes) - sum(outfluxes)
        return StateFlux(vcat(influxes, outfluxes), state, expr=state_expr)
    end

    function StateFlux(
        states::Pair{Num,Num};
    )
        ori_state, new_state = states[1], states[2]
        #* Construct the default calculation formula: new state variable minus old state variable
        state_expr = new_state - ori_state
        return StateFlux([new_state], ori_state, expr=state_expr)
    end

    function StateFlux(
        flux_names::Vector{Symbol},
        state_name::Symbol,
        param_names::Vector{Symbol};
        state_func::Function
    )
        #* Create variables by names
        fluxes = [first(@variables $nm) for nm in flux_names]
        state = first(@variables $state_name)
        params = [first(@parameters $nm) for nm in param_names]
        state_expr = state_func(fluxes, params)
        return StateFlux(fluxes, state, params, expr=state_expr)
    end

    function StateFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        state_name::Symbol;
    )
        influx_names, outflux_names = flux_names[1], flux_names[2]
        state_input_names = vcat(influx_names, outflux_names)
        #* Construct the default calculation function: sum of input variables minus sum of output variables
        state_func = (i, p) -> sum(i[1:length(influx_names)]) - sum(i[length(influx_names)+1:length(influx_names)+length(outflux_names)])
        return StateFlux(state_input_names, state_name, Symbol[], state_func=state_func)
    end

    function StateFlux(
        state_names::Pair{Symbol,Symbol}
    )
        ori_state_name, new_state_name = state_names[1], state_names[2]
        #* Create variables by names
        ori_state = first(@variables $ori_state_name)
        new_state = first(@variables $new_state_name)
        return StateFlux(ori_state => new_state)
    end
end

"""
(flux::AbstractStateFlux)(input::Vector, pas::ComponentVector; kwargs...)
(flux::AbstractStateFlux)(input::Array, pas::ComponentVector; kwargs...)
(flux::AbstractStateFlux)(input::Matrix, pas::ComponentVector; kwargs...)

Callable functions for AbstractStateFlux, this function is used for avoid not method match error

These methods are not supported for direct invocation. State fluxes should be used within a SimpleFlux or other composite flux structure.

# Arguments
- `input`: Input data (Array, Vector, or Matrix)
- `pas::ComponentVector`: Parameter vector
- `kwargs...`: Additional keyword arguments

# Returns
Throws an error indicating that State Flux cannot be run directly.

# Note
State fluxes represent internal model states and are not meant to be called independently. 
They should be incorporated into larger flux structures like SimpleFlux for proper execution within a hydrological model.

"""
(::AbstractStateFlux)(::Vector, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using SimpleFlux to run"
(::AbstractStateFlux)(::Matrix, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using SimpleFlux to run"
(::AbstractStateFlux)(::Array, ::ComponentVector, ::AbstractVector{Symbol}; kwargs...) = @error "State Flux cannot run directly, please using SimpleFlux to run"


(::AbstractRouteFlux)(::Vector, ::ComponentVector; kwargs...) = @error "Abstract RouteFlux is not support for single timepoint"
(route::AbstractRouteFlux)(input::Matrix, pas::ComponentVector; kwargs...) = @error "Must be implemented by subtype"
function (route::AbstractRouteFlux)(input::Array, pas::ComponentVector, ptypes::AbstractVector{Symbol}; kwargs...)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    pytype_params = [pas[:params][ptype] for ptype in ptypes]
    pytype_initstates = [pas[:initstates][ptype] for ptype in ptypes]

    sols = map(eachindex(ptypes)) do (idx)
        tmp_pas = ComponentVector(params=pytype_params[idx], initstates=pytype_initstates[idx])
        node_sols = reduce(hcat, route(input[:, idx, :], tmp_pas))
        node_sols
    end
    sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
    return permutedims(sol_arr, (1, 3, 2))
end

"""
    DiscRouteFlux{rtype} <: AbstractDiscRouteFlux

Represents a discrete routing flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `states::Vector{Num}`: A vector of state variables.
- `infos::NamedTuple`: Contains metadata about the flux, including input, output, parameter, and state names.

# Type Parameters
- `rtype`: A symbol representing the specific type of discrete routing flux.

# Constructors
    DiscRouteFlux(input::Num, params::Vector{Num}, states::Vector{Num}, routetype::Symbol)

# Description
`DiscRouteFlux` is a structure that encapsulates a discrete routing flux calculation in a hydrological model. 
It is designed to represent various types of discrete routing processes, with the specific type indicated 
by the `rtype` parameter.

The structure automatically generates output variable names based on the input, and organizes 
the information about inputs, outputs, parameters, and states into the `infos` field.

This structure is particularly useful for representing discrete routing processes in hydrological models, 
where water is transferred from one point to another in the system using discrete time steps or methods.
"""
struct DiscRouteFlux{rtype} <: AbstractDiscRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function DiscRouteFlux(
        input::Num,
        params::Vector{Num};
        routetype::Symbol,
        output::Union{Num,Nothing}=nothing,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_names = Symbolics.tosymbol.(params, escape=false)
        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end
        #* Setup the name information of the hydroroutement
        infos = (input=[input_name], output=[output_name], param=param_names, state=Symbol[])

        return new{routetype}(
            [input],
            [output],
            params,
            infos
        )
    end
end

"""
    ContRouteFlux{rtype} <: AbstractContRouteFlux

Represents a continuous routing flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `states::Vector{Num}`: A vector of state variables.
- `infos::NamedTuple`: Contains metadata about the flux, including input, output, parameter, and state names.

# Type Parameters
- `rtype`: A symbol representing the specific type of continuous routing flux.

# Constructors
    ContRouteFlux(input::Num, params::Vector{Num}, states::Vector{Num}, routetype::Symbol)

# Description
`ContRouteFlux` is a structure that encapsulates a continuous routing flux calculation in a hydrological model. 
It is designed to represent various types of continuous routing processes, with the specific type indicated 
by the `rtype` parameter.

The structure automatically generates output variable names based on the input, and organizes 
the information about inputs, outputs, parameters, and states into the `infos` field.

This structure is particularly useful for representing continuous routing processes in hydrological models, 
where water is transferred from one point to another in the system using continuous time steps.
"""
struct ContRouteFlux{rtype} <: AbstractContRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function ContRouteFlux(
        input::Num,
        params::Vector{Num};
        routetype::Symbol,
        output::Union{Num,Nothing}=nothing,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_names = Symbolics.tosymbol.(params, escape=false)
        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end
        #* Setup the name information of the hydroroutement
        infos = (input=[input_name], output=[output_name], param=param_names, state=Symbol[])

        return new{routetype}(
            [input],
            [output],
            params,
            infos
        )
    end
end

"""
    UnitHydroFlux{solvetype} <: AbstractRouteFlux

Represents a unit hydrograph flux model for routing water through a hydrological system.

# Fields
- `inputs::Vector{Num}`: A vector of input variables (Num).
- `outputs::Vector{Num}`: A vector of output variables (Num).
- `params::Vector{Num}`: A vector of parameter variables (Num).
- `uhfunc::Function`: The unit hydrograph function.
- `infos::NamedTuple`: A named tuple containing information about inputs, outputs, parameters, and states.

# Constructor
    UnitHydroFlux(
        input::Num,
        param::Num,
        uhfunc::Function;
        solvetype::Symbol=:unithydro1
    )

Constructs a UnitHydroFlux object.

# Arguments
- `input::Num`: The input variable.
- `param::Num`: The parameter variable.
- `uhfunc::Function`: The unit hydrograph function.
- `solvetype::Symbol`: The solver type (default is `:unithydro1`).

# Description
UnitHydroFlux represents a unit hydrograph flux model for routing water through a hydrological system. 
It uses a unit hydrograph function to transform input flows into routed output flows.

The structure supports different solving methods, specified by the `solvetype` parameter. 
Currently, it implements the `:unithydro1` solver type, which uses a discrete problem approach 
to calculate the routed flow.

This flux model is particularly useful in hydrological modeling for representing the 
temporal distribution of runoff as it moves through a watershed. It can account for the 
lag and attenuation of flow as water travels through the system.

The `uhfunc` field holds the unit hydrograph function, which typically takes a parameter 
(like time) and returns weights that describe how an input is distributed over time in the output.

When called, the UnitHydroFlux object applies the unit hydrograph to the input flow series, 
effectively convolving the input with the unit hydrograph to produce the routed output flow.

This structure is designed to be flexible and can be integrated into larger hydrological models 
to represent various routing processes in different parts of a water system.

"""
struct UnitHydroFlux{solvetype} <: AbstractUnitHydroFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "unit hydrology function"
    uhfunc::Function
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function UnitHydroFlux(
        input::Num,
        param::Num,
        uhfunc::Function;
        output::Union{Num,Nothing}=nothing,
        solvetype::Symbol=:unithydro1,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_name = Symbolics.tosymbol(param, escape=false)
        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end
        #* Setup the name information of the hydroroutement
        infos = (input=[input_name], output=[output_name], param=[param_name], state=Symbol[])

        return new{solvetype}(
            [input],
            [output],
            [param],
            uhfunc,
            infos
        )
    end
end

(::UnitHydroFlux)(::Vector, ::ComponentVector; kwargs...) = @error "UnitHydroFlux is not support for single timepoint"

#* 一种是用构建Discrete problem的方式求解
function (flux::UnitHydroFlux{:unithydro1})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    #* 首先将lagflux转换为discrete problem
    function lag_prob(u, p, t)
        u = circshift(u, -1)
        u[end] = 0.0
        input_vec[Int(t)] .* p[:weight] .+ u
    end

    uh_weight = flux.uhfunc(pas[:params][flux.infos[:param][1]])
    prob = DiscreteProblem(lag_prob, input_vec[1] .* uh_weight, (1, length(input_vec)), ComponentVector(weight=uh_weight))
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    reshape(Array(sol)[1, :], 1, length(input_vec))
end

#* 一种是用构建稀疏矩阵的方式求和
function (flux::UnitHydroFlux{:unithydro2})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    uh_weight = flux.uhfunc(pas[:params][flux.infos[:param][1]])
    uh_result = [-(i - 1) => uh_wi .* input_vec for (i, uh_wi) in enumerate(uh_weight)]
    uh_sparse_matrix = spdiagm(uh_result...)
    sum_route = sum(uh_sparse_matrix, dims=2)[1:end-length(uh_weight)+1]
    reshape(sum_route, 1, length(input_vec))
end

function (uh::AbstractUnitHydroFlux)(input::Array, pas::ComponentVector; ptypes::AbstractVector{Symbol}, kwargs...)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    pytype_params = [pas[:params][ptype] for ptype in ptypes]

    sols = map(eachindex(ptypes)) do (idx)
        tmp_pas = ComponentVector(params=pytype_params[idx])
        node_sols = reduce(hcat, uh(input[:, idx, :], tmp_pas))
        node_sols
    end
    sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
    return permutedims(sol_arr, (1, 3, 2))
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
- `infos::NamedTuple`: Contains metadata about the flux, including input, output, and neural network parameter names.
- `nnios::NamedTuple`: Contains information about the neural network's input and output structure.

# Description
`NeuralFlux` is a structure that encapsulates a neural network-based flux calculation in a hydrological model. 
It combines symbolic mathematics with neural networks to create a flexible and powerful representation of complex hydrological processes.

The structure automatically handles the integration of the neural network into the symbolic framework, 
allowing for seamless incorporation of machine learning models into traditional hydrological equations.

This structure is particularly useful for representing complex, non-linear relationships in hydrological systems 
where traditional equations may be insufficient or unknown.
"""
struct NeuralFlux <: AbstractNeuralFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of neural network names (Symbol) and its variables (Num)"
    nnparam::Symbolics.Arr
    "flux expressions to descripe the formula of the output variable"
    expr::Symbolics.Arr{Num,1}
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "neural network  information: keys contains: input, output, param"
    infos::NamedTuple
    "neural network build information: keys contains: input, output, param"
    nnios::NamedTuple

    function NeuralFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        chain::Lux.AbstractExplicitContainerLayer,
    )
        #* Get input and output variables
        input_vars, output_vars = fluxes[1], fluxes[2]
        #* Get the neural network name (neural flux param name) and object
        chain_name = chain.name
        #* Initialize model parameter type for model parameter dimension definition
        init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
        init_params_axes = getaxes(init_params)

        #* Define parameter variables according to initialization parameters: Define type as Vector{parameter length}
        chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(init_params))
        #* Use Symbolics.array_term to define the slow-building parameter variables:
        #* when the model is called, it is rebuilt into the ComponentVector type according to
        #* the axes of `init_params` and the Vector type of the parameter as the calculation parameter input
        lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, init_params_axes, size=size(chain_params))

        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(input_vars, escape=false)
        output_names = Symbolics.tosymbol.(output_vars, escape=false)

        #* Constructing neural network input and output variables
        nn_input_name = Symbol(chain_name, :_input)
        nn_output_name = Symbol(chain_name, :_output)

        #* The input and output of the model can only be of type Symbolics.Arr{Num, 1},
        #* so it cannot be defined based on input_vars and output_vars
        nn_input = first(@variables $(nn_input_name)[1:length(input_names)])
        nn_output = first(@variables $(nn_output_name)[1:length(output_names)])

        #* Constructing a calculation expression based on a neural network
        flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)

        nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, init_params_axes))

        #* neuralflux infos
        infos = (input=input_names, output=output_names, param=Symbol[], nn=[chain_name])
        new(
            input_vars,
            output_vars,
            chain_params,
            flux_expr,
            nn_func,
            infos,
            (input=nn_input, output=nn_output, paramlen=length(init_params)),
        )
    end

    function NeuralFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        chain::Lux.AbstractExplicitContainerLayer,
    )
        input_names, output_names = flux_names[1], flux_names[2]

        input_vars = [first(@variables $input_name) for input_name in input_names]
        output_vars = [first(@variables $output_name) for output_name in output_names]

        return NeuralFlux(input_vars => output_vars, chain)
    end
end

function (flux::AbstractNeuralFlux)(input::Vector, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][flux.infos[:nn][1]]
    flux.func(input, nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::Matrix, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][flux.infos[:nn][1]]
    flux.func(input', nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::Array, pas::ComponentVector, ::AbstractVector{Symbol})
    nn_params = pas[:nn][flux.infos[:nn][1]]
    #* array dims: (ts_len * node_names * var_names)
    flux_output_vec = [flux.func(input[:, i, :], nn_params) for i in 1:size(input)[2]]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (1, 3, 2))
end