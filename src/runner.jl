mutable struct HydroRunner{runtype} <: AbstractRunner
    model::AbstractComponent
    const_params::ComponentVector
    model_states::ComponentVector

    function HydroRunner(;
        model::AbstractComponent,
        runtype::Symbol,
        const_params::ComponentVector,
        model_states::ComponentVector,
    )
        return new{runtype}(
            model,
            const_params,
            model_states
        )
    end
end

"""
代表是一次性输出所有结果
"""
function (runner::HydroRunner{:single_run})(
    input::NamedTuple,
    params::ComponentVector,
)
    # parameters and initstates merge 

    # use the model

    # processes output
end

"""
代表是逐步输出计算结果
"""
function (runner::HydroRunner{:multi_run})(
    input::NamedTuple,
    params::ComponentVector,
)
    # parameters and initstates merge 
    a = input["a"]
    # use the model
    b = input["a"]

    # processes output
end

"""
代表是持续计算输出
"""
function (runner::HydroRunner{:conti_run})(
    input::NamedTuple,
    params::ComponentVector,
)
    # parameters and initstates merge 

    # use the model

    # processes output
end