struct EnsembleComponent{C1<:AbstractComponent,C2<:AbstractComponent} <: AbstractHydroWrapper
    ensemble_component::C1
    aggregate_component::C2
    ensemble_size::Integer
    meta::HydroMeta

    function EnsembleComponent(
        ensemble_component::C1,
        aggregate_component::C2,
        ensemble_size::Integer
    ) where {C1<:AbstractComponent,C2<:AbstractComponent}
        hybrid_meta = merge(ensemble_component.meta, aggregate_component.meta)

        new{C1,C2}(ensemble_component, aggregate_component, ensemble_size, ensemble_component.meta)
    end
end