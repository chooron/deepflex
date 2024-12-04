# Core Structs in HydroModels.jl

HydroModels.jl provides several core structs to represent different components of hydrological modeling. Each struct serves a specific purpose in the modeling process. Let's explore these core structs:

- Flux Component: designed to represent simple hydrological fluxes.
- Bucket Component: designed to represent simple hydrological bucket or lumped model.
- Route Component: designed to represent routing method in distributed or vector model.
- Model Component: designed to contain multiple components and connect them together and represent a complex model.

this components are subtype of the AbstractComponent struct. that this component contain the same infos, including the input, output, parameter, state and neuralnetwork infos.

moreover, to make the component more flexible, all of the component support the callable feature. which means you can obtain the output of components through function calls.

```julia
function (comp::AbstractComponent)(input::Vector, pas::ComponentVector; kwargs...)
    ...
end

function (comp::AbstractComponent)(input::Matrix, pas::ComponentVector; kwargs...)
    ...
end

function (comp::AbstractComponent)(input::Array, pas::ComponentVector; ptypes::AbstractVector{Symbol}, kwargs...)
    ...
end
```

the callable function of can be divide into three different input types:

- Vector input: vector input is used for single time step, the dimension of the input is (`num_variables`, ), this is support for ode problems solve.
- Matrix input: Matrix input is used for a single node input (lumped model), by input the matrix (`num_variables`, `time_steps`), the component can calculate the output for each time step of the ouput variables. this is support for lumped model. need to notice that the route component is not support matrix input (since the route component is designed for multiple nodes input).
- Array input: Array input is used for multiple nodes input (distributed model), by input the array (`num_variables`, `num_nodes`, `time_steps`), the component can calculate the output for each time step of the ouput variables for each node. this is support for distributed model.

all of the callable function should input the pas (parameter, init states, sometime also include the neuralnetwork params) together with the input data, the type is `ComponentVector` based on the `ComponentArrays.jl`. for example:

```julia
pas = ComponentVector(params = (...), initstates = (...), nn = (...))
```

while when input pas for multiple nodes (eg. distributed model), the pas should be change its struct as follows:

```julia
pas = ComponentVector(params = (ptype1 = (...), ptype2 = (...), ...), initstates = (ptype1 = (...), ptype2 = (...), ...), nn = (...))
```

where `ptype1, ptype2` are the node type in the distributed model, for example, in a catchment, we can have different type of land use/land cover (eg. forest, grassland, water body, wetland, etc.), each type of land cover can have different parameters, the pas for each node should include the parameters for all of the land cover type, and the parameters for each land cover type is stored in the `params` field of the pas.
while for `nn`, it can be shared by all of the nodes, which means the nn params is the same for all of the nodes.

for array input, the callable function need input a extra argument `ptypes::AbstractVector{Symbol}`, which is the node type in the distributed model, and the order of the node type should be the same as the order of the nodes in the input array.

for example, if we have 2 type nodes in the catchment, the input array should be (num_variables, 2, time_steps) and the ptypes should be `[:ptype1, :ptype2]`.
