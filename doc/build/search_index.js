var documenterSearchIndex = {"docs":
[{"location":"tutorials/run_a_bucket.html#HydroModels-Bucket-Model-Implementation-Example","page":"Run a Bucket Model","title":"HydroModels Bucket Model Implementation Example","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html#Overview","page":"Run a Bucket Model","title":"Overview","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"This document demonstrates the implementation and usage of a hydrological bucket model using the HydroModels framework. The code showcases both single-node and multi-node simulations using the ExpHydro model structure.","category":"page"},{"location":"tutorials/run_a_bucket.html#Dependencies","page":"Run a Bucket Model","title":"Dependencies","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"using CSV\nusing DataFrames\nusing ComponentArrays\nusing BenchmarkTools\nusing HydroModels\nusing ModelingToolkit","category":"page"},{"location":"tutorials/run_a_bucket.html#Model-Configuration","page":"Run a Bucket Model","title":"Model Configuration","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html#Parameter-Setup","page":"Run a Bucket Model","title":"Parameter Setup","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"The model uses the following parameters for the hydrological simulation:","category":"page"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"f: 0.01674478 (Infiltration parameter)\nSmax: 1709.461015 (Maximum soil water storage)\nQmax: 18.46996175 (Maximum discharge)\nDf: 2.674548848 (Degree-day factor)\nTmax: 0.175739196 (Maximum temperature threshold)\nTmin: -2.092959084 (Minimum temperature threshold)","category":"page"},{"location":"tutorials/run_a_bucket.html#Initial-States","page":"Run a Bucket Model","title":"Initial States","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"Initial conditions for the model:","category":"page"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"snowpack: 0.0\nsoilwater: 1303.004248","category":"page"},{"location":"tutorials/run_a_bucket.html#Data-Input","page":"Run a Bucket Model","title":"Data Input","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"The model uses time series data from a CSV file located at \"data/exphydro/01013500.csv\" containing:","category":"page"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"Day length (dayl)\nMean temperature (tmean)\nPrecipitation (prcp)","category":"page"},{"location":"tutorials/run_a_bucket.html#Implementation-Examples","page":"Run a Bucket Model","title":"Implementation Examples","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html#1.-Single-Node-Simulation","page":"Run a Bucket Model","title":"1. Single Node Simulation","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"# Setup input data\ninput = (lday=df[ts, \"dayl(day)\"], temp=df[ts, \"tmean(C)\"], prcp=df[ts, \"prcp(mm/day)\"])\nsolver = HydroModels.ManualSolver{true}()\nconfig = (solver=solver,)\n\n# Convert input to required format\ninput_arr = Matrix(reduce(hcat, collect(input[ele.meta.inputs]))')\n\n# Run simulation\nresults = ele(input_arr, pas, config=config, convert_to_ntp=true)","category":"page"},{"location":"tutorials/run_a_bucket.html#2.-Multi-Node-Simulation","page":"Run a Bucket Model","title":"2. Multi-Node Simulation","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"# Setup for multiple nodes\nnode_num = 10\nnode_names = [Symbol(:node, i) for i in 1:node_num]\n\n# Create parameter and state vectors for all nodes\nnode_params = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([params], length(node_names))))\nnode_initstates = ComponentVector(NamedTuple{Tuple(node_names)}(repeat([init_states], length(node_names))))\nnode_pas = ComponentVector(params=node_params, initstates=node_initstates)\n\n# Prepare input data for multiple nodes\ninput_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))\nnode_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))\nnode_input = permutedims(node_input, (2, 3, 1))\n\n# Run simulation with multiple nodes\nrun_kwgs = (ptypes=node_names, timeidx=ts)\nresult = ele(node_input, node_pas, kwargs=run_kwgs)","category":"page"},{"location":"tutorials/run_a_bucket.html#Model-Structure","page":"Run a Bucket Model","title":"Model Structure","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"The model uses a bucket structure defined in exphydro.jl with two main components:","category":"page"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"Surface water component (bucket_1):\nHandles precipitation partitioning (rainfall/snowfall)\nComputes potential evapotranspiration\nManages snowmelt processes\nSoil water component:\nManages soil water storage\nComputes actual evaporation\nGenerates baseflow and surface flow","category":"page"},{"location":"tutorials/run_a_bucket.html#Usage-Notes","page":"Run a Bucket Model","title":"Usage Notes","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"The code demonstrates flexibility in handling both single-node and multi-node simulations\nInput data should be properly formatted with required columns (dayl, tmean, prcp)\nParameters and initial states can be adjusted based on specific catchment characteristics\nThe model uses a manual solver for time-stepping","category":"page"},{"location":"tutorials/run_a_bucket.html#Time-Series-Processing","page":"Run a Bucket Model","title":"Time Series Processing","text":"","category":"section"},{"location":"tutorials/run_a_bucket.html","page":"Run a Bucket Model","title":"Run a Bucket Model","text":"The example processes 10,000 time steps (ts = collect(1:10000)) and can be adjusted based on data availability and simulation requirements.","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Run-A-ExpHydro-Model","page":"Run ExpHydro Model","title":"Run A ExpHydro Model","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html#Overview","page":"Run ExpHydro Model","title":"Overview","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"本文档将介绍如何使用已构建的exphdyro模型用于集总式(Single HRU)的计算和分布式(Multiple HURs)的计算","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Dependencies","page":"Run ExpHydro Model","title":"Dependencies","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"First, lets import the required packages.","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"using CSV\nusing DataFrames\nusing ComponentArrays\nusing BenchmarkTools\nusing NamedTupleTools\nusing Plots","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Model-Setup","page":"Run ExpHydro Model","title":"Model Setup","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"接着需要对模型的参数和初始状态进行设置","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Parameter-Configuration","page":"Run ExpHydro Model","title":"Parameter Configuration","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"Exphydro模型包含以下6个参数,我们需要使用ComponentVector来定义这些参数","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"Parameter Value Description\nf 0.01674478 Infiltration parameter\nSmax 1709.461015 Maximum soil water storage\nQmax 18.46996175 Maximum discharge\nDf 2.674548848 Degree-day factor\nTmax 0.175739196 Maximum temperature threshold\nTmin -2.092959084 Minimum temperature threshold","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"params = ComponentVector(f=0.01674478, Smax=1709.461015, Qmax=18.46996175,\n Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Initial-States","page":"Run ExpHydro Model","title":"Initial States","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"然后就是定义模型的初始状态,针对exphydro模型的两个计算模块分别设置其对应的状态变量snowpack和soilwater初始值","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"inistates = ComponentVector(snowpack=0.0, soilwater=1303.004248)","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Data-Preparation","page":"Run ExpHydro Model","title":"Data Preparation","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"The test uses hydrometeorological data from \"data/exphydro/01013500.csv\":","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"file_path = \"data/exphydro/01013500.csv\"\ndata = CSV.File(file_path)\ndf = DataFrame(data)\nts = collect(1:10000)  # Time series length","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"Input variables include:","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"Day length (dayl)\nMean temperature (tmean)\nPrecipitation (prcp)","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Model-Testing","page":"Run ExpHydro Model","title":"Model Testing","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html#Single-Node-Testing","page":"Run ExpHydro Model","title":"Single Node Testing","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"# Configure solver\nsolver = HydroModels.ManualSolver()\n\n# Run model with performance benchmarking\nresult = model(input, pas, config=(solver=solver, timeidx=ts), convert_to_ntp=true)\n\n# Visualize results\nplot(result.flow)\nplot!(df[ts, \"flow(mm)\"])","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Performance-Benchmarking-For-Single-Node-(by-using-BenchmarkTools.jl)","page":"Run ExpHydro Model","title":"Performance Benchmarking For Single Node (by using BenchmarkTools.jl)","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"@btime model(input, pas, config=(solver=solver, timeidx=ts), convert_to_ntp=true);","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Multi-Node-Testing-(Optional)","page":"Run ExpHydro Model","title":"Multi-Node Testing (Optional)","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"The code includes commented sections for multi-node testing:","category":"page"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"# Setup multiple nodes\nnode_num = 10\ninputs = repeat([input], node_num)\nptypes = [Symbol(:node, i) for i in 1:node_num]\n\n# Configure parameters and states for multiple nodes\nparams_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([params], node_num)))\ninit_states_multi = ComponentVector(NamedTuple{Tuple(ptypes)}(repeat([init_states], node_num)))\npas_multi = ComponentVector(params=params_multi, initstates=init_states_multi)\n\n# Run multi-node simulation\nresults = model(inputs, pas_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true)","category":"page"},{"location":"tutorials/run_a_exphydro_model.html#Performance-Benchmarking-For-Multi-Node-(by-using-BenchmarkTools.jl)","page":"Run ExpHydro Model","title":"Performance Benchmarking For Multi-Node (by using BenchmarkTools.jl)","text":"","category":"section"},{"location":"tutorials/run_a_exphydro_model.html","page":"Run ExpHydro Model","title":"Run ExpHydro Model","text":"@btime model(inputs, pas_multi, config=(solver=solver, timeidx=ts), convert_to_ntp=true);","category":"page"},{"location":"index.html#HydroModels.jl","page":"Home","title":"HydroModels.jl","text":"","category":"section"},{"location":"index.html#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"HydroModels.jl is a modern hydrological modeling framework that extends and enhances SUPERFLEX's design philosophy. Built on Julia language and the SciML (Scientific Machine Learning) ecosystem, it combines flexible model construction with computational efficiency, particularly supporting deep learning integration in hydrological modeling.","category":"page"},{"location":"index.html#Key-Features","page":"Home","title":"Key Features","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Flexible Model Construction: Supports development of lumped, semi-distributed, and distributed hydrological models\nDeep Learning Integration: Enables neural network integration for enhanced flux calculations and dynamic parameter estimation\nComputational Efficiency: Leverages Julia's high-performance capabilities and the SciML ecosystem\nGradient-Based Optimization: Supports advanced parameter optimization techniques\nComprehensive Framework: Provides tools for both traditional hydrological modeling and modern machine learning approaches","category":"page"},{"location":"index.html#Framework-Capabilities","page":"Home","title":"Framework Capabilities","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"HydroModels.jl offers:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Easy implementation and customization of hydrological models\nIntegration with Julia's scientific computing ecosystem\nSupport for various modeling approaches:\nTraditional conceptual models\nNeural network enhanced models\nDistributed hydrological systems","category":"page"},{"location":"index.html#Case-Studies","page":"Home","title":"Case Studies","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"The framework has been validated through various applications, including:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"M50 Model Implementation\nGR4J-based Distributed Hydrological Model","category":"page"},{"location":"index.html#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Check out our tutorials to start using HydroModels.jl:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Run a Bucket Model\nRun ExpHydro Model","category":"page"},{"location":"index.html#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"HydroModels\")","category":"page"},{"location":"index.html#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"HydroModels.jl is an open-source project available on Github. We welcome contributions from the community to help advance deep learning applications in hydrology.","category":"page"},{"location":"index.html#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"For detailed information about using HydroModels.jl, please refer to our comprehensive documentation and tutorials in this site.","category":"page"}]
}
