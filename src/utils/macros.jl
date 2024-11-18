using ModelingToolkit: isparameter, get_variables
using ModelingToolkit
using ModelingToolkit: Differential as D
using Symbolics
using BenchmarkTools
using Lux
using LuxCore
include("../HydroModels.jl")
HydroFlux = HydroModels.HydroFlux
NeuralFlux = HydroModels.NeuralFlux
StateFlux = HydroModels.StateFlux

macro hydroflux(eq::Expr)
    quote
        # Get left and right sides of equation
        lhs = $(eq.args[2])
        rhs = $(eq.args[3])

        # Extract output variables from left side
        outputs = lhs isa Vector ? lhs : [lhs]

        # Extract inputs and parameters from right side
        rhs_vars = get_variables(rhs)
        lhs_vars = Num.(get_variables(lhs))
        inputs = Num.(filter(x -> !isparameter(x), rhs_vars))
        params = Num.(filter(x -> isparameter(x), rhs_vars))

        # Build the HydroFlux constructor call

        HydroFlux(inputs, outputs, params, exprs=[$(eq.args[3])])
    end
end

function (chain::LuxCore.AbstractExplicitContainerLayer)(var::Vector{Num})
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
    nn_input = first(@variables $(nn_input_name)[1:length(inputs)])

    #* Constructing a calculation expression based on a neural network
    flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)

    return (chain=chain, input=var, expr=flux_expr, params=chain_params)
end

macro nnflux(eq::Expr)
    quote
        lhs = $(eq.args[2])
        chaininfo = $(eq.args[3])
        inputs = chaininfo.input isa Vector ? chaininfo.input : [chaininfo.input]
        outputs = lhs isa Vector ? lhs : [lhs]
        NeuralFlux(inputs, outputs, chaininfo.params, chaininfo.expr, chaininfo.chain)
    end
end

macro stateflux(eq::Expr)
    quote
        # Get left and right sides of equation
        lhs = $(eq.args[2])
        rhs = $(eq.args[3])

        # Extract output variables from left side
        outputs = lhs isa Vector ? lhs : [lhs]

        # Extract state variables from right side
        @assert length(outputs) == 1 "State flux should have only one state variable"

        # Extract inputs and parameters from right side
        rhs_vars = get_variables(rhs)
        lhs_vars = Num.(get_variables(lhs))
        inputs = Num.(filter(x -> !isparameter(x), rhs_vars))
        params = Num.(filter(x -> isparameter(x), rhs_vars))

        # Build the HydroFlux constructor call
        StateFlux(inputs, outputs[1], params, exprs=$(eq.args[3]))
    end
end