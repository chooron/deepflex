function Flow(input_names::Vector{Symbol}; parameters::ComponentVector=ComponentVector())
    SimpleFlux(
        input_names,
        [:Flow],
        parameters,
        flow_func
    )
end

function flow_func(
    input::gen_namedtuple_type([:Baseflow,:Surfaceflow], T),
    parameters::Nothing=nothing
)::gen_namedtuple_type([:Flow], T) where {T<:Number}
    (Flow=input[:Baseflow] .+ input[:Surfaceflow],)
end

function flow_func(
    input::gen_namedtuple_type([:Rainfall,:Percolation,:Saturation], T),
    parameters::Nothing=nothing
)::gen_namedtuple_type([:Flow], T) where {T<:Number}
    (Flow=input[:Rainfall] .+ input[:Percolation] .- input[:Saturation],)
end
