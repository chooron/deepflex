// #set text(font:("SimSun","Times New Roman"), fallback: false, size: 10pt)
#set align(start + top)

#align(center, text(17pt)[
  *HydroModels.jl: A Flexible and Efficient Framework for Hydrological Modeling*
])

= Abstract
this need to introduce the abstract content of the paper

#set heading(numbering: "1.1")
= Introduction
Hydrological modeling has undergone significant evolution over the past decades, reflecting our growing understanding of water systems and advancements in computational capabilities. This progression can be broadly categorized into three main phases: conceptual hydrological models, deep learning-based models, and the current hybrid approach combining deep learning with physical principles.

Conceptual hydrological models, developed in the mid-20th century, aimed to represent watershed processes using simplified mathematical equations. These models, while effective in many scenarios, often struggled with complex, non-linear hydrological processes and required extensive calibration. As computational power increased, distributed and semi-distributed models emerged, allowing for more detailed spatial representation of watershed characteristics and processes.

The advent of machine learning, particularly deep learning, in the early 21st century brought a paradigm shift to hydrological modeling. Deep learning models, leveraging large datasets and powerful algorithms, demonstrated remarkable capabilities in capturing complex, non-linear relationships in hydrological systems. These data-driven approaches often outperformed traditional models in prediction tasks but lacked the physical interpretability crucial for understanding underlying processes.

The current frontier in hydrological modeling lies in hybrid approaches that combine the strengths of deep learning with physical principles. These physics-informed neural networks and theory-guided data science methods aim to leverage the predictive power of machine learning while maintaining consistency with known physical laws and domain knowledge. This approach promises more robust, interpretable, and generalizable models capable of handling the complexities of hydrological systems across diverse conditions.

The construction of hydrological models has evolved alongside these conceptual advancements. Conceptual models typically involve defining interconnected components representing various hydrological processes, often requiring careful parameterization and calibration. Distributed and semi-distributed models necessitate more complex structures to account for spatial heterogeneity, often incorporating GIS data and advanced routing schemes. Deep learning models, in contrast, rely on designing appropriate network architectures and training procedures to learn from large datasets.

Despite these advancements, the development of hydrological models remains challenging. The diversity of modeling approaches has led to the creation of various model building frameworks, each with its own strengths and limitations. For instance, the Framework for Understanding Structural Errors (FUSE) provides a flexible approach to conceptual hydrological modeling by allowing the combination of different model structures. The Structure for Unifying Multiple Modeling Alternatives (SUMMA) offers a unified approach to hydrological modeling, enabling the representation of multiple modeling approaches within a single framework. Superflex provides a highly flexible framework for conceptual hydrological modeling, allowing for the creation of customized model structures. The Modular Assessment of Rainfall–Runoff Models Toolbox (MARRMoT) offers a collection of conceptual hydrological models and tools for model comparison and evaluation. However, these frameworks often focus on specific types of models or modeling approaches, and may lack the flexibility to easily incorporate machine learning techniques or hybrid approaches. Moreover, the complexity of modern hybrid approaches can make model construction and modification a daunting task, requiring expertise in both hydrology and machine learning, which is not always fully addressed by existing frameworks.

These challenges highlight the need for a unified, flexible framework for hydrological model construction. Such a framework should accommodate various modeling paradigms, from simple conceptual models to complex hybrid approaches, while providing a consistent interface for model development, calibration, and analysis. It should facilitate the integration of different modeling components, allow for easy experimentation with novel approaches, and promote reproducibility in hydrological research.

In response to these needs, we present HydroModels.jl, a comprehensive and adaptable framework designed to streamline the process of hydrological model development and application. By providing a unified platform for diverse modeling approaches, HydroModels.jl aims to accelerate innovation in hydrological science and enhance our ability to address critical water resource challenges in an era of environmental change.

= Architecture of HydroModels.jl

HydroModels.jl is built upon a modular and flexible architecture, designed to accommodate a wide range of hydrological modeling approaches. The framework's core structure comprises four main classes: Flux, Bucket, Route, and Model. Each of these classes plays a crucial role in representing different aspects of the hydrological system, from individual processes to complete model structures.

== Flux Class

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

== Bucket Class

The Bucket class is a fundamental component in HydroModels.jl, representing water storage modules within the hydrological system. It serves as a versatile abstraction for various water reservoirs, including soil moisture, groundwater, and surface water bodies, enabling the simulation of dynamic changes in water storage over time.

At its core, the Bucket class dynamically generates functions for solving ordinary differential equations (ODEs) and processing non-StateFlux components. These functions are constructed at runtime using anonymous functions, efficiently integrating multiple flux components. This approach allows for flexible and efficient representation of hydrological processes within the bucket system. The detailed implementation of this process will be discussed later in the paper.

The reason for using this approach is that this method is highly efficient and quick, avoiding matrix update calculations and other computational overheads. It allows for flexible and efficient representation of the relationships between different hydrological processes within the bucket, enabling seamless integration of various flux components and adaptability to diverse modeling scenarios.

The Bucket class seamlessly integrates multiple flux components, manages state variables representing water storage, and implements ODEs for dynamic simulations. This versatility allows modelers to customize the bucket to represent various types of water storage units, adapting to the specific needs of different hydrological scenarios. The primary implementation, HydroBucket, exemplifies this flexibility, offering a customizable model that can be tailored to a wide range of water storage representations.

== Route Class

The Route class plays a crucial role in HydroModels.jl by simulating water movement through landscapes and river networks. It serves as the key differentiator between lumped, multi-node, semi-distributed, and fully distributed hydrological models. While these models may share similar runoff generation processes, their routing calculations vary significantly, leading to distinct model types.

The Route class is extended into three main subclasses, each catering to different spatial representations and routing approaches:

1. SumRoute: This subclass is designed for multi-node models, integrating results from multiple calculation units. It implements a simple weighted accumulation routing scheme, allowing for the aggregation of outputs from various model components.

2. GridRoute: Tailored for fully distributed models, GridRoute specializes in routing processes for grid-based calculation units. It primarily operates on flow direction matrices of watersheds, enabling detailed spatial representation of water movement across a gridded landscape.

3. VectorRoute: This subclass is particularly suited for semi-distributed models, focusing on routing processes in watershed network calculations. It employs directed graph computations to represent and simulate water flow through river networks and channel systems.

The functionality of these Route subclasses, particularly GridRoute and VectorRoute, is enhanced by their unique ability to support different RouteFlux implementations. This feature sets HydroModels.jl apart from many other modeling frameworks. Both GridRoute and VectorRoute construct a unified system of ordinary differential equations that incorporates two types of states: those specific to the RouteFlux components and those pertaining to the Route itself. This approach allows for a more comprehensive and flexible representation of routing processes, enabling the integration of various hydrological phenomena such as time delays and flow attenuation. The detailed implementation of this unified ODE system is provided in the appendix.

The Route class and its subclasses offer several advantages in hydrological modeling. They enable the construction of models with varying levels of spatial complexity, from simple lumped models to sophisticated distributed systems. This flexibility allows researchers and practitioners to choose the most appropriate spatial representation for their specific modeling needs, balancing computational efficiency with the desired level of detail in representing hydrological processes.

== Model Class

The Model class represents the highest-level structure in HydroModels.jl, serving as an abstract base for various hydrological system representations. Its primary implementation, HydroModel, integrates multiple Flux, Bucket, and Route components to create comprehensive watershed models. HydroModel capitalizes on the uniform interfaces of these underlying components, enabling a streamlined and efficient simulation process.

HydroModel orchestrates the overall simulation by sequentially iterating through its constituent components. This approach leverages the standardized interfaces of Flux, Bucket, and Route classes, The consistency of these interfaces allows for flexible combination of modules within the framework while maintaining conceptual integrity. This design philosophy supports a wide range of model configurations, enabling researchers to construct and experiment with diverse hydrological representations tailored to specific research needs or watershed characteristics, all within a coherent and unified modeling environment. Further, the sequential computation simplifies model execution and facilitates easy integration of diverse hydrological processes and spatial representations.

By managing the data flow between components, handling time stepping, and generating outputs, the Model class provides a unified framework for hydrological modeling. This design allows HydroModels.jl to accommodate a wide range of modeling approaches, from simple lumped models to complex, spatially-distributed systems. The architecture's emphasis on modularity and extensibility enables researchers and practitioners to easily implement new processes or modeling techniques within the existing structure, fostering innovation and adaptability in hydrological modeling.

The Model class thus serves as a powerful tool for developing, testing, and applying diverse hydrological models. Its ability to seamlessly integrate various components while maintaining a consistent computational approach underscores the flexibility and robustness of HydroModels.jl as a comprehensive hydrological modeling framework.

= Methodologies and Fundamental Algorithms of HydroModels.jl

== Construction Methods

The construction methods are designed to transform user-defined components and parameters into cohesive, computationally efficient structures that represent complex hydrological systems.

=== Symbolic Programming and Runtime Function Generation

The construction methods in HydroModels.jl are fundamentally based on symbolic programming, metadata extraction, and runtime function generation. This approach forms the backbone of the framework's flexibility and efficiency.

At the core of HydroModels.jl's construction methods is symbolic programming, implemented using the Symbolics.jl package, which is analogous to Python's SymPy. Symbolic programming serves as the foundation for building various types within the framework, providing crucial indicators for inputs, outputs, and parameters of different components.

In the context of hydrological modeling, symbolic programming allows for the representation of intermediate fluxes in the model's computational process. These symbolic variables, when combined with the framework's core classes, enable the representation of complex hydrological processes. Essentially, these variables act as intermediate fluxes in the model, while the classes define how these fluxes are calculated within the model structure.

Following the symbolic representation, the framework generates metadata for each component. This metadata, stored in string format, records essential information about each class, including its inputs, outputs, parameters, and states. This feature allows users to access component information readily and facilitates efficient data computation. The metadata serves as a guide for the framework, enabling quick access to component characteristics and aiding in the orchestration of data flow within the model.

The final step in the construction process involves the building of anonymous functions for specific components. This process varies depending on the type of component:

1. Flux and Bucket Components: For runoff generation modules, anonymous functions are typically constructed at runtime. At its core, the Bucket class dynamically generates two types of functions: one for solving ordinary differential equations (ODEs) and another for processing non-StateFlux components. These functions are constructed at runtime using anonymous functions, which efficiently sequence and integrate multiple flux components. The functions are built once during the model initialization phase and remain fixed thereafter, avoiding repeated construction during model execution. The function generation process utilizes the `build_ele_func` method, leveraging the `RuntimeGeneratedFunctions.jl` and `SymbolicUtils.jl` packages. This approach employs symbolic computation techniques to seamlessly combine various flux components, including neural network-based fluxes, into a cohesive and efficient computational structure. The resulting functions process input data to output all internal variables of the fluxes, enabling a flexible and powerful representation of hydrological processes within the bucket system.

2. RouteFlux and Route Components: The construction of these components depends on specific requirements. RouteFlux components are first constructed based on the chosen routing method (e.g., Muskingum routing, Nash unit hydrograph). Then, the Route component is built and its routing calculation method is set. These components do not typically require runtime anonymous function construction.

3. UnitHydroFlux Components: For these components, the main requirement is to specify the type of unit hydrograph (e.g., triangular unit hydrograph) to be used in the calculations.

This construction approach offers several advantages, including simplicity in building, strong developmental potential, and high computational efficiency, making HydroModels.jl a powerful and flexible framework for hydrological modeling.

== Computational Methods

HydroModels.jl employs advanced computational strategies to enhance model parameter calibration and computational efficiency, particularly for multi-unit simulations. The framework supports both single-node (lumped models) and multi-node (distributed models) computations through innovative array-based approaches. The key computational methods are as follows:

=== Dimensional Handling for Input Arrays

In model simulations, input arrays are categorized based on the number of calculation units:

1. Two-dimensional matrices (N × T): Represent single calculation unit inputs, where N is the variable dimension and T is the time length.
2. Three-dimensional arrays (N × M × T): Represent multi-unit inputs, where M is the number of calculation units.

=== Metadata-Driven Data Extraction and Result Storage

HydroModels.jl employs a sophisticated metadata-driven approach for efficient data extraction and result storage during model execution. This method enhances computational efficiency and facilitates seamless data flow between model components. The key aspects of this approach are:

1. Metadata Generation and Pre-computation Planning: During model construction, metadata for each component is generated and stored, including information about inputs, outputs, and intermediate fluxes. Based on this metadata, the framework pre-determines the sequence of calculations, optimizing the computational flow before simulations begin.

2. Efficient Data Extraction and Result Storage: During execution, the framework uses the pre-stored metadata to extract precisely the required input data for each module and store results in a predetermined order. This approach minimizes unnecessary data processing and allows for efficient concatenation of results.

3. Optimized Data Flow: By leveraging the metadata-driven approach, the framework ensures smooth data flow between components without redundant calculations or data reorganization, managing intermediate fluxes based on the pre-planned sequence.

This metadata-driven method significantly enhances the framework's efficiency, particularly in complex models with multiple interconnected components, allowing for streamlined data management and improved overall model performance.

=== Slicing and Broadcasting Techniques

Since the generated computation functions support only single time point calculations, the framework implements array broadcasting and slicing techniques for both 2D and 3D array computations:

1. For 2D matrices (single calculation unit):
   - Time dimension slicing creates sub-arrays (an iterator).
   - Parameters are broadcast to match the dimensions of the sliced iterator.
   - Results are computed for each time step using the computation function.

2. For 3D arrays (multiple calculation units):
   - Slicing occurs along both time and calculation unit dimensions.
   - Both input arrays and parameters are sliced for each calculation unit.
   - Computation proceeds similarly to the 2D matrix case for each unit.

The array broadcasting mechanism automatically matches input parameters with sliced arrays, enabling element-wise computation and producing output results of appropriate dimensions.

=== ODE Solving for Multiple Calculation Units

To ensure simultaneous solving across all calculation units in ODE computations, the framework represents equation state variables as 2D matrices. Each matrix represents the state variables for all calculation units at a specific time point. This approach allows for efficient parallel solving of ODEs across multiple units.

By supporting both 2D and 3D array computations, HydroModels.jl provides a flexible and efficient framework for a wide range of hydrological modeling applications, from simple lumped models to complex distributed systems. This computational approach enhances the framework's ability to handle diverse spatial representations and computational requirements in hydrological modeling.

== Optimization and Solving Methods

HydroModels.jl leverages the powerful SciML ecosystem, particularly the DifferentialEquations and Optimization packages, to provide robust optimization and solving capabilities. This deep integration allows users to access a wide range of state-of-the-art methods while maintaining flexibility and extensibility.

=== Solving Methods

The framework supports both continuous and discrete ordinary differential equations (ODEs), corresponding to ODEProblem and DiscreteProblem in the DifferentialEquations library. Users can select from a variety of solvers and adjust precision settings to suit their specific modeling needs. This flexibility allows for efficient handling of diverse hydrological processes, from rapid surface runoff to slow groundwater movements.

For continuous ODEs, users can choose from methods such as Runge-Kutta, Adams, or BDF (Backward Differentiation Formula) solvers, depending on the stiffness and characteristics of their system. Discrete ODEs, often used in conceptual models, can be solved using appropriate discrete solvers.

=== Optimization Approaches

HydroModels.jl offers two main optimization approaches: black-box optimization and gradient-based optimization. Both methods construct objective functions that measure the discrepancy between model outputs and observed data.

1. Black-box Optimization: This approach treats the model as a black box, making it suitable for complex models or when derivatives are difficult to compute. It utilizes algorithms that do not require gradient information, such as Particle Swarm Optimization or Differential Evolution.

2. Gradient-based Optimization: Leveraging the Zygote automatic differentiation library and SciMLSensitivity, this method computes gradients of the model computations and ODE solutions. It allows for more efficient optimization, especially for large-scale problems. However, it requires careful implementation to ensure compatibility with Zygote's computation rules throughout the entire calculation process.

The choice between these optimization methods depends on the specific model characteristics, computational resources, and the desired balance between optimization speed and robustness. By providing these diverse optimization capabilities, HydroModels.jl enables researchers to effectively calibrate their models and explore parameter spaces, enhancing the overall model performance and reliability.

== Features of HydroModels.jl

=== Modularity and Composability:
The framework is built on a highly modular architecture, allowing for the seamless integration of various hydrological components. The primary building blocks include:

- Flux components (SimpleFlux, StateFlux, NeuralFlux, RouteFlux)
- Bucket models (HydroBucket)
- Routing schemes (WeightSumRoute, GridRoute, VectorRoute)
- Comprehensive models (HydroModel)

This modular design enables researchers and practitioners to construct complex hydrological models by combining these components in various configurations, tailoring the model structure to specific research needs or watershed characteristics.

=== Abstraction and Extensibility:
HydroModels.jl employs a system of abstract types (e.g., AbstractComponent, AbstractFlux, AbstractRoute) to provide a clear hierarchy and enable easy extensibility. This design allows users to implement new components or modify existing ones without altering the core framework, promoting adaptability to diverse hydrological scenarios.

=== Integration of Traditional and Machine Learning Approaches:
The framework seamlessly incorporates both traditional hydrological equations and machine learning techniques. The NeuralFlux component, for instance, allows for the integration of neural networks within the hydrological modeling process, enabling hybrid modeling approaches that can capture complex, non-linear relationships in hydrological systems.

=== Efficient Data Flow and Computation:
HydroModels.jl is designed with performance in mind. The framework utilizes efficient data structures and algorithms to manage the flow of information between components. For example, the HydroModel struct uses input indices to efficiently map overall model inputs to component-specific inputs, optimizing the simulation process.

=== Support for Multiple Spatial Representations:
The framework accommodates various spatial representations of hydrological systems, from lumped models to distributed grid-based and vector-based approaches. This is evident in the different routing schemes provided (WeightSumRoute, GridRoute, VectorRoute), allowing for flexible spatial modeling of water movement.

=== Advanced Numerical Methods:
HydroModels.jl incorporates sophisticated numerical methods for solving differential equations and optimizing parameters. The framework includes utilities for parameter optimization and supports various ODE solvers, enabling accurate and efficient simulation of hydrological processes over time.

=== Metadata Management:
The framework places a strong emphasis on metadata management. Each component maintains detailed information about its inputs, outputs, states, and parameters, facilitating model introspection, documentation, and debugging.

=== Flexibility in Model Application:
HydroModels.jl supports both single-node and multi-node simulations, allowing for applications ranging from simple catchment models to complex, spatially distributed hydrological systems. The framework's design enables easy scaling from local to regional modeling efforts.

=== Interoperability:
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
