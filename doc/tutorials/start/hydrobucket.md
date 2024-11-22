# Introduction to HydroBucket Structure

HydroBucket is a fundamental component in hydrological modeling, representing a computational unit of a hydrological module. It consists of simple fluxes and state fluxes.

As an independent hydrological module, HydroBucket can be used to represent either a lumped hydrological model or a specific computational module within a larger hydrological model (such as the snowmelt calculation module in the Exp-hydro model).

## Core function of the HydroBucket

### 1. runtime function generate

`HydroBucket` integrates multiple simple fluxes and state fluxes to describe flux calculations and state balance calculations for a hydrological module. Although individual flux types can produce outputs based on data inputs, repeatedly calling the callable functions of fluxes within HydroBucket presents several challenges:

1. Pre-allocating a result matrix: When flux calculations replace values at corresponding positions, and these values are directly accessed when needed, this approach, while computationally efficient, involves mutating arrays. This violates the constraints of the Zygote.jl optimizer used in model optimization (see limitations of Zygote.jl).

2. To avoid mutating arrays, we initially used NamedTuples to store flux calculation results during development, with each output name corresponding to a flux calculation result. Although less efficient than the previous method, this approach avoids mutating arrays and satisfies the Zygote.jl optimizer constraints for model optimization.

3. To enhance computational efficiency while meeting the Zygote.jl optimizer constraints, the HydroBucket constructor uses the expressions provided in simple fluxes and state fluxes to generate a runtime function (based on SymbolicUtils.jl) during model construction. This function efficiently computes all flux results (flux_func) and supports the construction and rapid solution of ODE problems (ode_func).

**The detailed implementation of runtime function generation**

Runtime function generation primarily relies on the Symbolics.jl and SymbolicUtils.jl packages. The function generates a runtime function based on the variables, parameters, and expressions stored in the fluxes. The specific steps are as follows:

1. Establish relationships between output variables and their respective expressions using the Assignment function from SymbolicUtils.jl.

2. Organize output variables into an Array using MakeArray.

3. Construct function parameters using `[DestructuredArgs(inputs), DestructuredArgs(params), DestructuredArgs(nns)]`. This generates runtime functions with input parameter formats:
   - `function fluxes(input::Vector, params::Vector, nns::Vector)`
   - `function dfluxes(input::Vector, state::Vector, params::Vector, nns::Vector)`

4. Build the function body content using Let, which includes both the relationships between output variables and their expressions, and the organization of output variables into an Array.

5. Generate the function expression using Func, and create the runtime function using the `@RuntimeGenerated` macro.


### 2. ODE problem solve

When the input argument `dfluxes` is not empty, HydroBucket constructs the `ode_func` runtime function and uses the `OrdinaryDiffEq.jl` package to build and solve an ODE problem. The process involves the following steps:

1. Linear interpolation is used to construct interpolation functions for input fluxes. This allows for obtaining flux values at any time point within the calculation period when solving the ODE problem.

2. To accommodate the requirements of Optimization.jl, we use ComponentVector for parameter input. We pre-extract the index describing the parameter order, enabling quick indexing of corresponding parameter values during flux calculations, thus improving computational efficiency.

3. Finally, we have built wrapper structs based on the `solve` function from the `OrdinaryDiffEq.jl` package. These wrappers `Solver` are used for both continuous and discrete ODE solving.