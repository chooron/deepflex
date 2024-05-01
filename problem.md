I'm using Optimization.jl to optimize some parameters of ODEProblem, where AutoForwardDiff is used for adtype, when the code is executed for the first time, the program does not show any errors, but when I execute the program again without restarting julia, I get the following error:

```
ERROR: MethodError: no method matching Float64(::ForwardDiff.Dual{ForwardDiff.Tag{OptimizationForwardDiffExt.var"#37#55"{…}, Float64}, Float64, 6})

Closest candidates are:
  (::Type{T})(::Real, ::RoundingMode) where T<:AbstractFloat
   @ Base rounding.jl:207
  (::Type{T})(::T) where T<:Number
   @ Core boot.jl:792
  Float64(::UInt8)
   @ Base float.jl:165
  ...

Stacktrace:
  [1] convert
    @ D:\software\Julia-1.10.0\pkg\packages\ForwardDiff\PcZ48\src\dual.jl:433 [inlined]
  [2] setindex!(A::Vector{ForwardDiff.Dual{…}}, x::ForwardDiff.Dual{ForwardDiff.Tag{…}, ForwardDiff.Dual{…}, 6}, i1::Int64)
    @ Base .\array.jl:1021
  [3] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\SymbolicUtils\c0xQb\src\code.jl:418 [inlined]
  [4] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\Symbolics\Eas9m\src\build_function.jl:546 [inlined]
  [5] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\SymbolicUtils\c0xQb\src\code.jl:375 [inlined]
  [6] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\RuntimeGeneratedFunctions\M9ZX8\src\RuntimeGeneratedFunctions.jl:163 [inlined]
  [7] macro expansion
    @ .\none:0 [inlined]
  [8] generated_callfunc
    @ .\none:0 [inlined]
  [9] (::RuntimeGeneratedFunctions.RuntimeGeneratedFunction{…})(::Vector{…}, ::Vector{…}, ::Vector{…}, ::Float64)
    @ RuntimeGeneratedFunctions D:\software\Julia-1.10.0\pkg\packages\RuntimeGeneratedFunctions\M9ZX8\src\RuntimeGeneratedFunctions.jl:150
 [10] (::ModelingToolkit.var"#f#709"{…})(du::Vector{…}, u::Vector{…}, p::ModelingToolkit.MTKParameters{…}, t::Float64)
    @ ModelingToolkit D:\software\Julia-1.10.0\pkg\packages\ModelingToolkit\kByuD\src\systems\diffeqs\abstractodesystem.jl:344
 [11] (::SciMLBase.Void{ModelingToolkit.var"#f#709"{…}})(::Vector{ForwardDiff.Dual{…}}, ::Vararg{Any})
    @ SciMLBase D:\software\Julia-1.10.0\pkg\packages\SciMLBase\QEvkv\src\utils.jl:482
 [12] (::FunctionWrappers.CallWrapper{…})(f::SciMLBase.Void{…}, arg1::Vector{…}, arg2::Vector{…}, arg3::ModelingToolkit.MTKParameters{…}, arg4::Float64)
    @ FunctionWrappers D:\software\Julia-1.10.0\pkg\packages\FunctionWrappers\Q5cBx\src\FunctionWrappers.jl:65
 [13] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\FunctionWrappers\Q5cBx\src\FunctionWrappers.jl:137 [inlined]
 [14] do_ccall
    @ D:\software\Julia-1.10.0\pkg\packages\FunctionWrappers\Q5cBx\src\FunctionWrappers.jl:125 [inlined]
 [15] FunctionWrapper
    @ D:\software\Julia-1.10.0\pkg\packages\FunctionWrappers\Q5cBx\src\FunctionWrappers.jl:144 [inlined]
 [16] _call
    @ D:\software\Julia-1.10.0\pkg\packages\FunctionWrappersWrappers\9XR0m\src\FunctionWrappersWrappers.jl:12 [inlined]
 [17] FunctionWrappersWrapper
    @ D:\software\Julia-1.10.0\pkg\packages\FunctionWrappersWrappers\9XR0m\src\FunctionWrappersWrappers.jl:10 [inlined]
 [18] ODEFunction
    @ D:\software\Julia-1.10.0\pkg\packages\SciMLBase\QEvkv\src\scimlfunctions.jl:2296 [inlined]
 [19] initialize!(integrator::OrdinaryDiffEq.ODEIntegrator{…}, cache::OrdinaryDiffEq.Rosenbrock23Cache{…})
    @ OrdinaryDiffEq D:\software\Julia-1.10.0\pkg\packages\OrdinaryDiffEq\ZbQoo\src\perform_step\rosenbrock_perform_step.jl:10
 [20] __init(prob::ODEProblem{…}, alg::OrdinaryDiffEq.Rosenbrock23{…}, timeseries_init::Tuple{}, ts_init::Tuple{}, ks_init::Tuple{}, recompile::Type{…}; saveat::Float64, tstops::Tuple{}, d_discontinuities::Tuple{}, save_idxs::Nothing, save_everystep::Bool, save_on::Bool, save_start::Bool, save_end::Nothing, callback::Nothing, dense::Bool, calck::Bool, dt::Float64, dtmin::Float64, dtmax::Float64, force_dtmin::Bool, adaptive::Bool, gamma::Rational{…}, abstol::Float64, reltol::Float64, qmin::Rational{…}, qmax::Int64, qsteady_min::Int64, qsteady_max::Rational{…}, beta1::Nothing, beta2::Nothing, qoldinit::Rational{…}, controller::Nothing, fullnormalize::Bool, failfactor::Int64, maxiters::Int64, internalnorm::typeof(DiffEqBase.ODE_DEFAULT_NORM), internalopnorm::typeof(LinearAlgebra.opnorm), isoutofdomain::typeof(DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN), unstable_check::typeof(DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK), verbose::Bool, timeseries_errors::Bool, dense_errors::Bool, advance_to_tstop::Bool, stop_at_next_tstop::Bool, initialize_save::Bool, progress::Bool, progress_steps::Int64, progress_name::String, progress_message::typeof(DiffEqBase.ODE_DEFAULT_PROG_MESSAGE), progress_id::Symbol, userdata::Nothing, allow_extrapolation::Bool, initialize_integrator::Bool, alias_u0::Bool, alias_du0::Bool, initializealg::OrdinaryDiffEq.DefaultInit, kwargs::@Kwargs{})
    @ OrdinaryDiffEq D:\software\Julia-1.10.0\pkg\packages\OrdinaryDiffEq\ZbQoo\src\solve.jl:518
 [21] __init (repeats 5 times)
    @ D:\software\Julia-1.10.0\pkg\packages\OrdinaryDiffEq\ZbQoo\src\solve.jl:11 [inlined]
 [22] #__solve#761
    @ D:\software\Julia-1.10.0\pkg\packages\OrdinaryDiffEq\ZbQoo\src\solve.jl:6 [inlined]
 [23] __solve
    @ D:\software\Julia-1.10.0\pkg\packages\OrdinaryDiffEq\ZbQoo\src\solve.jl:1 [inlined]
 [24] solve_call(_prob::ODEProblem{…}, args::OrdinaryDiffEq.Rosenbrock23{…}; merge_callbacks::Bool, kwargshandle::Nothing, kwargs::@Kwargs{…})
    @ DiffEqBase D:\software\Julia-1.10.0\pkg\packages\DiffEqBase\NaUtB\src\solve.jl:612
 [25] solve_up(prob::ODEProblem{…}, sensealg::SciMLSensitivity.ForwardSensitivity{…}, u0::Vector{…}, p::ModelingToolkit.MTKParameters{…}, args::OrdinaryDiffEq.Rosenbrock23{…}; kwargs::@Kwargs{…})
    @ DiffEqBase D:\software\Julia-1.10.0\pkg\packages\DiffEqBase\NaUtB\src\solve.jl:1080
 [26] solve_up
    @ D:\software\Julia-1.10.0\pkg\packages\DiffEqBase\NaUtB\src\solve.jl:1066 [inlined]
 [27] solve(prob::ODEProblem{…}, args::OrdinaryDiffEq.Rosenbrock23{…}; sensealg::SciMLSensitivity.ForwardSensitivity{…}, u0::Nothing, p::Nothing, wrap::Val{…}, kwargs::@Kwargs{…})
    @ DiffEqBase D:\software\Julia-1.10.0\pkg\packages\DiffEqBase\NaUtB\src\solve.jl:1003
 [28] (::Main.DeepFlex.ODESolver)(ode_prob::ODEProblem{…}, state_names::Vector{…})
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\utils\solver.jl:13
 [29] (::Main.DeepFlex.HydroElement{…})(input::@NamedTuple{…}, pas::ComponentVector{…}; solver::Main.DeepFlex.ODESolver)
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\element.jl:90
 [30] (::Vector{…})(input::@NamedTuple{…}, pas::ComponentVector{…}; solver::Main.DeepFlex.ODESolver)
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\element.jl:148
 [31] (::Main.DeepFlex.HydroNode)(input::@NamedTuple{…}, pas::ComponentVector{…}; solver::Main.DeepFlex.ODESolver)
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\node.jl:54
 [32] (::Main.DeepFlex.HydroNode)(input::@NamedTuple{…}, pas::ComponentVector{…})
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\node.jl:46
 [33] (::Main.DeepFlex.var"#predict_func#39"{…})(x::Vector{…}, p::SciMLBase.NullParameters)
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\utils\optimize.jl:74
 [34] (::Main.DeepFlex.var"#objective#40"{…})(x::Vector{…}, p::SciMLBase.NullParameters)
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\utils\optimize.jl:77
 [35] (::OptimizationForwardDiffExt.var"#37#55"{…})(::Vector{…})
    @ OptimizationForwardDiffExt D:\software\Julia-1.10.0\pkg\packages\OptimizationBase\rRpJs\ext\OptimizationForwardDiffExt.jl:98
 [36] #39
    @ D:\software\Julia-1.10.0\pkg\packages\OptimizationBase\rRpJs\ext\OptimizationForwardDiffExt.jl:102 [inlined]
 [37] vector_mode_dual_eval!
    @ D:\software\Julia-1.10.0\pkg\packages\ForwardDiff\PcZ48\src\apiutils.jl:24 [inlined]
 [38] vector_mode_gradient!(result::Vector{…}, f::OptimizationForwardDiffExt.var"#39#57"{…}, x::Vector{…}, cfg::ForwardDiff.GradientConfig{…})
    @ ForwardDiff D:\software\Julia-1.10.0\pkg\packages\ForwardDiff\PcZ48\src\gradient.jl:96
 [39] gradient!
    @ ForwardDiff D:\software\Julia-1.10.0\pkg\packages\ForwardDiff\PcZ48\src\gradient.jl:37 [inlined]
 [40] (::OptimizationForwardDiffExt.var"#38#56"{…})(::Vector{…}, ::Vector{…})
    @ OptimizationForwardDiffExt D:\software\Julia-1.10.0\pkg\packages\OptimizationBase\rRpJs\ext\OptimizationForwardDiffExt.jl:102
 [41] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\OptimizationOptimisers\AOkbT\src\OptimizationOptimisers.jl:68 [inlined]
 [42] macro expansion
    @ D:\software\Julia-1.10.0\pkg\packages\Optimization\5DEdF\src\utils.jl:32 [inlined]
 [43] __solve(cache::OptimizationCache{…})
    @ OptimizationOptimisers D:\software\Julia-1.10.0\pkg\packages\OptimizationOptimisers\AOkbT\src\OptimizationOptimisers.jl:66
 [44] solve!(cache::OptimizationCache{…})
    @ SciMLBase D:\software\Julia-1.10.0\pkg\packages\SciMLBase\QEvkv\src\solve.jl:188
 [45] solve(::OptimizationProblem{…}, ::Adam; kwargs::@Kwargs{…})
    @ SciMLBase D:\software\Julia-1.10.0\pkg\packages\SciMLBase\QEvkv\src\solve.jl:96
 [46] solve
    @ D:\software\Julia-1.10.0\pkg\packages\SciMLBase\QEvkv\src\solve.jl:93 [inlined]
 [47] param_grad_optim(component::Main.DeepFlex.HydroNode; tunable_pas::ComponentVector{…}, const_pas::ComponentVector{…}, input::@NamedTuple{…}, target::@NamedTuple{…}, kwargs::@Kwargs{})
    @ Main.DeepFlex f:\julia\DeepFlex.jl\src\utils\optimize.jl:82
 [48] top-level scope
    @ f:\julia\DeepFlex.jl\test\optimization\test_node_grad_optimization.jl:38
Some type information was truncated. Use `show(err)` to see complete types.
```

This problem comes from a project I am currently writing, I have the above problem when I use ModelingToolkit.jl to build ODEProblem and use AutoForwardDiff as the Adtype, and I don't have this problem when I build ODEProblem directly through the Function or use AutoFiniteDiff, see the detailed execution file https://github.com/chooron/DeepFlex.jl/blob/main/test/optimization/test_node_grad_optimization.jl
