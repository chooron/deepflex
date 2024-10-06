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

function (runner::HydroRunner{:single_run})(
    input::NamedTuple,
    params::ComponentVector,
)
    # parameters and initstates merge 

    # use the model

    # processes output
end

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

function (runner::HydroRunner{:conti_run})(
    input::NamedTuple,
    params::ComponentVector,
)
    # parameters and initstates merge 

    # use the model

    # processes output
end