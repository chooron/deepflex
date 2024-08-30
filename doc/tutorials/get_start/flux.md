# Introduction to Flux Structure

Flux is a fundamental concept in hydrological modeling ,which denote a kind of hydrological process, such soil infiltration, surface runoff, and groundwater flow. In this pakage, we use flux to describe those fluxes in one hydrological model. based on different computation process, we proposed different types of flux:

## 1. Simple Flux

Simple flux is used to describe the process of common hydrological process that can be descripe by a simple computation formula, such as potential evaporation process:

\[
\text{pet} = 29.8 \times \text{lday} \times 24 \times 0.611 \times \exp\left(\frac{17.3 \times \text{temp}}{\text{temp} + 237.3}\right) / (\text{temp} + 273.2)
\]

In the potential evaporation formula, we can set `lday` and `temp` to represent the day length and temperature, and `pet` to represent the potential evaporation, and the `pet` can be calculated by the formula by inputting the `lday` and `temp`. to build this flux in `HydroModel`, we need to:

```julia
pet = SimpleFlux([lday, temp] => [pet], exprs = [lday * 24 * 0.611 * exp(17.3 * temp / (temp + 237.3)) / (temp + 273.2)])
```

In the code, `[lday, temp]` is the input variables, `[pet]` is the output variable, and `exprs` is the formula of the flux.

Then, we present an another example of simple flux that contain parameters:

In exp-hydro model, the melt can be calculated by the formula:

\[
\text{melt} = \text{step\_func}(\text{temp} - \text{Tmax}) \times \text{step\_func}(\text{snowpack}) \times \min(\text{snowpack}, \text{Df} \times (\text{temp} - \text{Tmax}))
\]

where `temp` is the temperature, `Tmax` is the maximum temperature, `Df` is the degree-day factor, `snowpack` is the snowpack, and `melt` is the melt, and the `step_func` is the step function, which is used to limit the melt when the temperature is greater than the maximum temperature or the snowpack is zero.

To build this flux in `HydroModel`, we need to:

```julia
melt = SimpleFlux([temp, Tmax, snowpack, Df] => [melt], [Tmax, Df],
    exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))])
```

In the code, `[temp, Tmax, snowpack, Df]` is the input variables, `[melt]` is the output variable, `[Tmax, Df]` is the parameters, and `exprs` is the formula of the flux.

As we can observe when constructing a SimpleFlux, both the output variables and expressions are of array type. This indicates that we can use a single SimpleFlux instance to represent multiple formulas. However, it's crucial to ensure that the variables required in the calculation formulas correspond to the input variables.

For example, if we use `a` and `b` to calculate `c`, but the calculation of `d` requires `c` as an input, we would need to create two separate SimpleFlux entities. This separation ensures that the dependencies between variables are properly maintained and that each flux represents a distinct step in the calculation process.

Through multiple SimpleFlux instances, we can easily represent a hydrological calculation module (excluding state variables). Of course, due to the flexible nature of HydroModels.jl, SimpleFlux can be directly input into the Model without the need to compose an additional bucket. This flexibility allows for efficient model construction and customization. For more details on how to incorporate SimpleFlux directly into models, please refer to the `Model` introduction section.

## Special simple flux: neural network flux (NeuralFlux)

`NeuralFlux` is a specialized structure in HydroModels.jl that integrates neural networks into hydrological flux calculations. For example, we can use a neural network to approximate the potential evaporation process:

\[
\text{pet} = \text{petnn}(\text{lday}, \text{temp})
\]

the implementation of this flux is as follows:

```julia
nn = Lux.Chain(
    Lux.Dense(2, 10),
    Lux.relu,
    Lux.Dense(10, 1),
    name=:petnn
)
pet = NeuralFlux([lday, temp] => [pet], nn)
```
In this code, we use `Lux.jl` to build a neural network (specified name by `name=:petnn`, it was used for the keys of the neural network parameters), and the neural network has two inputs (`lday` and `temp`) and one output (`pet`). It is important to note that the name of the input and output flux must be the same as the input dims and output dims in the neural network.

## 2. State Flux

State Flux represents a special and crucial component in hydrological models (or modules). This type of flux is based on the law of mass conservation, which states that the difference between the input and output fluxes of a module equals the change in its storage.

For example, the Exp-Hydro model contains two calculation modules. In the snowmelt module, the rate of change in snow depth is represented by the formula:

\[
\frac{d\text{snowpack}}{dt} = \text{snowfall} - \text{melt}
\]

The input flux of the module is snowfall, and the output flux is snowmelt. The change in storage is represented by the snow depth. This can be expressed using `StateFlux` as follows:

```julia
state_flux = StateFlux([snowfall] => [melt], snowpack)
```

In this example, the input flux is `snowfall`, the output flux is `melt`, and the storage is `snowpack`.

The construction method above is derived from a simple balance calculation based on pre-calculated fluxes. However, for some models where state variables are updated after flux calculations, or where the dstate calculation formula is non-linear, a more general construction method is needed:

```julia
state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr)]
```

In this code, the calculation of \(\frac{dsoilwater}{dt}\) requires the use of functions like `exp` and `step_func` (some state calculations may even involve parameters). To accommodate this, we can adopt a construction approach similar to SimpleFlux. This means providing input variables, output variables (in this case, the state value), parameters, and expressions to implement more complex state calculation formulas.

This approach allows for greater flexibility in defining state dynamics, enabling the incorporation of non-linear functions and parameter-dependent calculations. By structuring StateFlux in this way, we can handle a wider range of hydrological processes and model complexities, making it easier to represent sophisticated state variable dynamics within the HydroModel.jl framework.

need to notice that, in the StateFlux, one state flux can only be used to describe one state variable, if you want to describe multiple state variables, you need to use multiple StateFlux.

## 3. Route Flux

Route Flux can be used for both routing process (kinematic wave, discharge model, muskingum method, etc.) in distributed hydrological model and routing function (like unit hydrograph, linear reservoir, etc.) in conceptual model. the calculation progress of this kind of flux is usually complex, which can not be described by a simple formula, and as for different routing methods, the calculation process is different, we categorize them into three types:

- `UnitHydroFlux`: This flux is specifically designed for the unit hydrograph method, which is commonly used in conceptual hydrological models. The calculation process of this flux is based on unit hydrograph theory, where we need to specify a unit-hydrograph function to generate the unit hydrograph weights. It's important to note that UnitHydroFlux calculations require input of the entire flux data and perform calculations uniformly. This means we cannot extract data for a single time point to output to other grid points, making it unsuitable for routing calculations in grid and vector route. The UnitHydroFlux is primarily used for lumped conceptual models or weighted sum route (we will introduce this in the route section).

- `GridRouteFlux`: This flux is specifically designed for grid routing, which is commonly used in grided hydrological models. The internal calculation process of this routing method is constructed as a continuous ordinary differential equation. It assumes that the flux flows continuously from one node to another, meaning the flux calculation process is continuous. For more details on the grid routing process, please refer to the [route](route.md) section.

- `VectorRouteFlux`: This flux is specifically designed for vector routing, which is commonly used in vectorized hydrological models. The internal calculation process of this routing method is constructed as a discrete calculation process. This means that the flux movement occurs at each calculation time step. If we were to construct an ODE (Ordinary Differential Equation) calculation process, it would need to be built for each time step, which is clearly impractical. Therefore, we adopt a discrete approach for construction. The flux moves from one node to another at each discrete time step, allowing for efficient computation in vector-based models. This approach is particularly suitable for routing methods like the unit hydrograph. For more details on the vector routing  process, please refer to the [route](route.md) section.


We have implemented some route functions in the `fluxes` directory, including discharge model, muskingum, linear reserivoir, and unit hydrograph. the information table of these fluxes is as follows:

| Flux Name | Description | Parameters | Routetype |
|-----------|-------------|------------|-----------|
| DischargeRouteFlux | Hydrological discharge model | lag, s_river | GridRouteFlux |
| CascadeRouteFlux | Nash Cascade (Linear Reservoir) | k, n | GridRouteFlux |
| MuskingumRouteFlux | Muskingum Method | k, x | VectorRouteFlux |
| UnitHydroFlux | Unit Hydrograph | lag_func | UnitHydroFlux |
