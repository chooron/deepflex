@with_kw_noshow mutable struct Unit <: Component
    elements::Vector{Element}
    topology::Dict
end