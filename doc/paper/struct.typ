#set text(font:("SimSun","Times New Roman"), fallback: false, size: 12pt)
#set align(start + top)

#align(center, text(17pt)[
  *A paper outline of HydroModels.jl*
])

= Abstract
this need to introduce the abstract content of the paper

#set heading(numbering: "1.1")
= Introduction
本内容需要介绍水文模型的发展,包括概念性水文模型到深度学习水文模型再到现在的深度学习+物理的水文模型;
然后需要介绍水文模型的搭建方法,包括概念式水文模型,(半)分布式水文模型,以及深度学习水文模型;
然后引出当前水文模型搭建的单一性,以及搭建的困难性,所以需要一个通用的水文模型搭建框架;

= Design Philosophy of HydroModels.jl

== Core Architecture of HydroModels.jl

HydroModels.jl is built upon a modular and flexible architecture, designed to accommodate a wide range of hydrological modeling approaches. The framework's core structure comprises four main classes: Flux, Bucket, Route, and Model. Each of these classes plays a crucial role in representing different aspects of the hydrological system, from individual processes to complete model structures.

=== Flux Class

The Flux class serves as the fundamental building block of HydroModels.jl, representing basic hydrological processes. This class embodies the mathematical formulations that govern the movement and transformation of water within various components of the hydrological cycle. By design, the Flux class accommodates a diverse array of water fluxes, ranging from precipitation and evapotranspiration to infiltration and runoff. Its strength lies in its ability to not only encapsulate process-specific equations but also efficiently manage inputs, outputs, and parameters. Moreover, the Flux class demonstrates remarkable versatility in supporting both straightforward algebraic relationships and more complex differential equations. This flexibility enables it to capture a broad spectrum of hydrological processes, spanning from simple linear interactions to intricate dynamic systems, thereby providing a robust framework for comprehensive hydrological modeling.

To accommodate various modeling needs, the Flux class is extended into several specialized subclasses, see the following table:

#table(
  columns: (auto, auto, auto),
  inset: 10pt,
  align: horizon,
  [*Flux Type*], [*Description*], [*Supported Dimensions*],
  [SimpleFlux], [Implements basic algebraic relationships for straightforward process representations.], [Vector, Matrix, Array],
  [StateFlux], [Designed for processes involving state variables, allowing for more complex dynamics.], [Vector, Matrix, Array],
  [NeuralFlux], [Incorporates machine learning techniques, enabling the integration of data-driven approaches within the traditional modeling framework.], [Vector, Matrix, Array],
  [RouteFlux], [Specialized for routing processes, handling water movement through the system.], [Vector, Matrix, Array],
  [UnitHydroFlux], [Implements the unit hydrograph concept for rainfall-runoff modeling.], [Vector, Matrix, Array],
  [TimeVaryingFlux], [Addresses time-dependent processes, allowing for time series inputs.], [Vector, Matrix, Array]
)

These diverse Flux subclasses, while tailored to specific hydrological processes, share a common interface and structure within the HydroModels.jl framework. This uniformity is a key strength of the design, ensuring consistency in metadata management, function calls, and overall usage across different flux types. Each subclass, regardless of its specialization, adheres to a standardized approach for handling inputs, outputs, and parameters. This consistency not only simplifies the implementation of new flux types but also enhances the framework's usability, allowing modelers to seamlessly integrate and interchange various flux components within their hydrological models. The unified structure facilitates a modular approach to model construction, where different flux types can be easily combined or substituted without necessitating significant changes to the overall model architecture.

=== Bucket Class

The Bucket class is a fundamental component in HydroModels.jl, representing water storage modules within the hydrological system. It serves as a versatile abstraction for various water reservoirs, including soil moisture, groundwater, and surface water bodies, enabling the simulation of dynamic changes in water storage over time.

At its core, the Bucket class dynamically generates two types of functions: one for solving ordinary differential equations (ODEs) and another for processing non-StateFlux components. These functions are constructed at runtime using anonymous functions, which efficiently sequence and integrate multiple flux components. The functions are built once during the model initialization phase and remain fixed thereafter, avoiding repeated construction during model execution. The function generation process utilizes the `build_ele_func` method, leveraging the `RuntimeGeneratedFunctions.jl` and `SymbolicUtils.jl` packages. This approach employs symbolic computation techniques to seamlessly combine various flux components, including neural network-based fluxes, into a cohesive and efficient computational structure. The resulting functions process input data to output all internal variables of the fluxes, enabling a flexible and powerful representation of hydrological processes within the bucket system. 

The reason for using this approach is that this method is highly efficient and quick, avoiding matrix update calculations and other computational overheads. It allows for flexible and efficient representation of the relationships between different hydrological processes within the bucket, enabling seamless integration of various flux components and adaptability to diverse modeling scenarios.

The Bucket class seamlessly integrates multiple flux components, manages state variables representing water storage, and implements ODEs for dynamic simulations. This versatility allows modelers to customize the bucket to represent various types of water storage units, adapting to the specific needs of different hydrological scenarios. The primary implementation, HydroBucket, exemplifies this flexibility, offering a customizable model that can be tailored to a wide range of water storage representations.

=== Route Class

The Route class plays a crucial role in HydroModels.jl by simulating water movement through landscapes and river networks. It serves as the key differentiator between lumped, multi-node, semi-distributed, and fully distributed hydrological models. While these models may share similar runoff generation processes, their routing calculations vary significantly, leading to distinct model types.

The Route class is extended into three main subclasses, each catering to different spatial representations and routing approaches:

1. SumRoute: This subclass is designed for multi-node models, integrating results from multiple calculation units. It implements a simple weighted accumulation routing scheme, allowing for the aggregation of outputs from various model components.

2. GridRoute: Tailored for fully distributed models, GridRoute specializes in routing processes for grid-based calculation units. It primarily operates on flow direction matrices of watersheds, enabling detailed spatial representation of water movement across a gridded landscape.

3. VectorRoute: This subclass is particularly suited for semi-distributed models, focusing on routing processes in watershed network calculations. It employs directed graph computations to represent and simulate water flow through river networks and channel systems.

The functionality of these Route subclasses, particularly GridRoute and VectorRoute, is enhanced by their unique ability to support different RouteFlux implementations. This feature sets HydroModels.jl apart from many other modeling frameworks. Both GridRoute and VectorRoute construct a unified system of ordinary differential equations that incorporates two types of states: those specific to the RouteFlux components and those pertaining to the Route itself. This approach allows for a more comprehensive and flexible representation of routing processes, enabling the integration of various hydrological phenomena such as time delays and flow attenuation. The detailed implementation of this unified ODE system is provided in the appendix.

The Route class and its subclasses offer several advantages in hydrological modeling. They enable the construction of models with varying levels of spatial complexity, from simple lumped models to sophisticated distributed systems. This flexibility allows researchers and practitioners to choose the most appropriate spatial representation for their specific modeling needs, balancing computational efficiency with the desired level of detail in representing hydrological processes.

=== Model Class

The Model class serves as the highest-level structure in HydroModels.jl, integrating multiple Flux, Bucket, and Route components to create complete hydrological models. It manages the overall simulation process, including data flow between components, time stepping, and output generation.

Key features of the Model class include:
- Composition of multiple hydrological components
- Management of input-output relationships between components
- Execution of the overall simulation process

The primary implementation of the Model class is:

- HydroModel: The main model class that integrates all components and orchestrates the simulation.

This modular design philosophy allows HydroModels.jl to accommodate a wide spectrum of hydrological modeling approaches, from simple lumped models to complex, spatially-distributed systems. The framework's architecture emphasizes extensibility, enabling users to easily implement new processes or modeling techniques within the existing structure. By providing these core classes and their specialized subclasses, HydroModels.jl offers a flexible and powerful tool for hydrological researchers and practitioners to develop, test, and apply a diverse range of hydrological models.


== Features of HydroModels.jl

== Modularity and Composability:
The framework is built on a highly modular architecture, allowing for the seamless integration of various hydrological components. The primary building blocks include:

- Flux components (SimpleFlux, StateFlux, NeuralFlux, RouteFlux)
- Bucket models (HydroBucket)
- Routing schemes (WeightSumRoute, GridRoute, VectorRoute)
- Comprehensive models (HydroModel)

This modular design enables researchers and practitioners to construct complex hydrological models by combining these components in various configurations, tailoring the model structure to specific research needs or watershed characteristics.

== Abstraction and Extensibility:
HydroModels.jl employs a system of abstract types (e.g., AbstractComponent, AbstractFlux, AbstractRoute) to provide a clear hierarchy and enable easy extensibility. This design allows users to implement new components or modify existing ones without altering the core framework, promoting adaptability to diverse hydrological scenarios.

== Integration of Traditional and Machine Learning Approaches:
The framework seamlessly incorporates both traditional hydrological equations and machine learning techniques. The NeuralFlux component, for instance, allows for the integration of neural networks within the hydrological modeling process, enabling hybrid modeling approaches that can capture complex, non-linear relationships in hydrological systems.

== Efficient Data Flow and Computation:
HydroModels.jl is designed with performance in mind. The framework utilizes efficient data structures and algorithms to manage the flow of information between components. For example, the HydroModel struct uses input indices to efficiently map overall model inputs to component-specific inputs, optimizing the simulation process.

== Support for Multiple Spatial Representations:
The framework accommodates various spatial representations of hydrological systems, from lumped models to distributed grid-based and vector-based approaches. This is evident in the different routing schemes provided (WeightSumRoute, GridRoute, VectorRoute), allowing for flexible spatial modeling of water movement.

== Advanced Numerical Methods:
HydroModels.jl incorporates sophisticated numerical methods for solving differential equations and optimizing parameters. The framework includes utilities for parameter optimization and supports various ODE solvers, enabling accurate and efficient simulation of hydrological processes over time.

== Metadata Management:
The framework places a strong emphasis on metadata management. Each component maintains detailed information about its inputs, outputs, states, and parameters, facilitating model introspection, documentation, and debugging.

== Flexibility in Model Application:
HydroModels.jl supports both single-node and multi-node simulations, allowing for applications ranging from simple catchment models to complex, spatially distributed hydrological systems. The framework's design enables easy scaling from local to regional modeling efforts.

== Interoperability:
The framework is designed to be interoperable with other Julia packages and external tools. It leverages Julia's ecosystem for tasks such as data interpolation, graph computations, and deep learning, enhancing its capabilities and ease of use.

In conclusion, the design philosophy of HydroModels.jl emphasizes flexibility, modularity, and efficiency, providing a powerful tool for hydrological modeling. By combining traditional hydrological concepts with modern computational techniques, HydroModels.jl offers a versatile platform for researchers and practitioners to develop, test, and apply a wide range of hydrological models, from simple conceptual representations to complex, spatially-distributed systems.

== Appendix

=== build_ele_func

1. Variable Preparation:
   - Collects variables from all input functions (funcs and dfuncs).
   - Converts these variables to symbols and creates a named tuple.

2. Function Parameter Construction:
   - Extracts input, state, and parameter variables from the prepared variables.
   - Collects neural network parameters (if present).

3. Assignment and Output List Construction:
   - Iterates through each input function, handling both regular and neural fluxes.
   - For neural fluxes, it matches inputs to neural network inputs and outputs to calculation results.
   - For regular fluxes, it matches output variables to their corresponding expressions.
   - Builds a comprehensive list of assignments and outputs.

4. Function Generation:
   - Constructs argument lists for both flux and state functions.
   - Uses `@RuntimeGeneratedFunction` to create efficient, dynamically generated functions.
   - For flux functions, it combines all flux outputs into a single vector function.
   - For state functions (if present), it creates a separate function for differential equations.

5. Output:
   - Returns two functions: a merged flux function and a merged state function (if applicable).
