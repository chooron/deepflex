using ModelingToolkit: isparameter, get_variables
using Symbolics

macro hydroflux(eq::Expr)
    quote
        # Get left and right sides of equation
        lhs = $(eq.args[2])
        rhs = $(eq.args[3])
        println(lhs)
        println(rhs)

        # Extract output variables from left side
        outputs = lhs isa Vector ? lhs : [lhs]

        # Extract inputs and parameters from right side
        rhs_vars = get_variables(rhs)
        lhs_vars = Num.(get_variables(lhs))
        inputs = Num.(filter(x -> !isparameter(x), rhs_vars))
        params = Num.(filter(x -> isparameter(x), rhs_vars))

        # Build the HydroFlux constructor call

        NamedTuple(
            [:inputs => inputs,
            :outputs => lhs_vars,
            :params => params]
        )
    end
end


@variables a b c d
@parameters k1 k2
@hydroflux a ~ k1 * b + k2 * c + d
