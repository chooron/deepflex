"""
$(TYPEDEF)
A struct representing common hydrological fluxes
# Fields
$(FIELDS)
# Example
```
# define a common function
function flow_func(
    i::namedtuple(:baseflow, :surfaceflow),
    p::NamedTuple;
    kw...
)
    i[:baseflow] .+ i[:surfaceflow]
end

flow_flux = SimpleFlux(
    [:baseflow, :surfaceflow],
    :flow,
    param_names=Symbol[],
    func=flow_func
)
```
"""
struct SimpleFlux <: AbstractSimpleFlux
    "A map of input names (Symbol) and its variables (Num)"
    input_info::NamedTuple
    "A map of output names (Symbol) and its variables (Num)"
    output_info::NamedTuple
    "A map of parameters names (Symbol) and its variables (Num)"
    param_info::NamedTuple
    "flux expressions to descripe the formula of the output variable"
    flux_exprs::Vector{Num}

    function SimpleFlux(
        input_info::NamedTuple,
        output_info::NamedTuple,
        param_info::NamedTuple,
        flux_exprs::Vector{Num}
    )
        return new(
            input_info,
            output_info,
            param_info,
            flux_exprs
        )
    end

    function SimpleFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        param_names::Vector{Symbol}=Symbol[];
        flux_funcs::Vector{<:Function}=Function[],
        kwargs...
    )
        #* Get input and output names
        input_names, output_names = flux_names[1], flux_names[2]

        if length(flux_funcs) > 0
            #* Create variables by names
            inputs = [first(@variables $var) for var in input_names]
            outputs = [first(@variables $var) for var in output_names]
            params = [first(@parameters $var) for var in param_names]
            #* When a calculation function is provided, exprs are constructed based on the calculation function and variables
            flux_exprs = [flux_func(inputs, params) for flux_func in flux_funcs]
        else
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            inputs, outputs, params = hydro_equation.inputs, hydro_equation.outputs, hydro_equation.params
            flux_exprs = HydroEquations.expr(hydro_equation; kwargs...)
        end

        #* Building the struct
        return SimpleFlux(
            NamedTuple{Tuple(input_names)}(inputs),
            NamedTuple{Tuple(output_names)}(outputs),
            NamedTuple{Tuple(param_names)}(params),
            flux_exprs
        )
    end

    function SimpleFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        flux_exprs::Vector{Num}
    )
        #* Get input and output variables
        inputs, outputs = fluxes[1], fluxes[2]

        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = Symbolics.tosymbol.(params)

        return SimpleFlux(
            NamedTuple{(Tuple(input_names))}(inputs),
            NamedTuple{(Tuple(output_names))}(outputs),
            NamedTuple{(Tuple(param_names))}(params),
            flux_exprs,
        )
    end
end

struct StateFlux <: AbstractStateFlux
    "A map of input names (Symbol) and its variables (Num)"
    input_info::NamedTuple
    "A map of output names (Symbol) and its variables (Num)"
    output_info::NamedTuple
    "A map of parameter names (Symbol) and its variables (Num)"
    param_info::NamedTuple
    "flux expressions to descripe the formula of the state variable"
    state_expr::Num

    function StateFlux(
        fluxes::Vector{Num},
        state::Num,
        params::Vector{Num};
        state_expr::Num
    )
        #* Convert to a symbol based on the variable
        state_input_names = Symbolics.tosymbol.(fluxes, escape=false)
        state_name = Symbolics.tosymbol(state, escape=false)
        state_param_names = Symbolics.tosymbol.(params, escape=false)
        return new(
            NamedTuple{Tuple(state_input_names)}(fluxes),
            NamedTuple{tuple(state_name)}([state]),
            NamedTuple{Tuple(state_param_names)}(params),
            state_expr,
        )
    end

    function StateFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        state::Num;
    )
        influxes, outfluxes = fluxes[1], fluxes[2]
        #* Construct the default calculation formula: sum of input variables minus sum of output variables
        state_expr = sum(influxes) - sum(outfluxes)
        return StateFlux(vcat(influxes, outfluxes), state, Num[], state_expr=state_expr)
    end

    function StateFlux(
        states::Pair{Num,Num};
    )
        ori_state, new_state = states[1], states[2]
        #* Construct the default calculation formula: new state variable minus old state variable
        state_expr = new_state - ori_state
        return StateFlux([new_state], ori_state, Num[], state_expr=state_expr)
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
        return StateFlux(fluxes, state, params, state_expr=state_expr)
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
$(TYPEDEF)
A hydrological flux calculated via a neural network (based on `Lux.jl`)
# Fields
$(FIELDS)
# Example
```
et_ann = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)
etnn_flux = NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, param_names=:etnn, chain=et_ann)
```
"""
struct NeuralFlux <: AbstractNeuralFlux
    "A map of input names (Symbol) and its variables (Num)"
    input_info::NamedTuple
    "A map of output names (Symbol) and its variables (Num)"
    output_info::NamedTuple
    "nn input and output information"
    nn_info::NamedTuple
    "flux expression"
    flux_expr

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

        new(
            NamedTuple{Tuple(input_names)}(input_vars),
            NamedTuple{Tuple(output_names)}(output_vars),
            (name=chain_name, input=nn_input, output=nn_output, params=chain_params),
            flux_expr
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

"""
$(TYPEDEF)
A flux used in hydrological unit-hydrograph calculations
# Fields
$(FIELDS)
# Example
```
slowflow_lagflux = LagFlux(:slowflow, :slowflow_lag, lag_func=LumpedHydro.uh_1_half, lag_time=:x4)
```
"""
struct LagFlux <: AbstractLagFlux
    "A map of input names (Symbol) and its variables (Num)"
    input_info::NamedTuple
    "A map of output names (Symbol) and its variables (Num)"
    output_info::NamedTuple
    "A map of parameter names (Symbol) and its variables (Num)"
    param_info::NamedTuple
    "lag functions"
    lag_func::Function
    
    function LagFlux(
        flux_names::Pair{Symbol,Symbol},
        lag_time_name::Symbol,
        lag_func::Function;
        kwargs...,
    )
        #* Create variables by names
        fluxes = [first(@variables $nm) for nm in flux_names]
        lag_time = first(@parameters $lag_time_name)

        new(
            NamedTuple{(flux_names[1],)}([fluxes[1]]),
            NamedTuple{(flux_names[2],)}([fluxes[2]]),
            NamedTuple{(lag_time_name,)}([lag_time]),
            lag_func
        )
    end

    function LagFlux(
        fluxes::Pair{Num,Num},
        lag_time::Num,
        lag_func::Function;
        kwargs...,
    )
        #* Convert to a symbol based on the variable
        flux_name_1 = Symbolics.tosymbol(fluxes[1], escape=false)
        flux_name_2 = Symbolics.tosymbol(fluxes[2], escape=false)
        lag_time_name = Symbolics.tosymbol(lag_time, escape=false)

        new(
            NamedTuple{(flux_name_1,)}([fluxes[1]]),
            NamedTuple{(flux_name_2,)}([fluxes[2]]),
            NamedTuple{(lag_time_name,)}([lag_time]),
            lag_func
        )
    end
end