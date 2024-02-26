"""
RoutingStore in GR4J
"""
function RoutingStore_GR4J(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:Q9]),
        Baseflow([:RoutingStore], parameters=parameters[[:x3, :γ]]),
        Recharge([:RoutingStore], parameters=parameters[[:x2, :x3, :ω]]),
        Flow([:Baseflow, :Recharge, :Q1])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(RoutingStore=input[:Q9] + input[:Recharge] - input[:Baseflow])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

function RoutingStore_HYMOD(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Splitter([:Saturation], [:Pf, :Ps], parameters=ComponentVector(Pf=parameters[:a], Ps=1 - parameters[:a])),
        SimpleFlux([:F1], [:Qf1], parameters[[:kf]], (i, p) -> p[:kf] * i[:F1]),
        SimpleFlux([:F2], [:Qf2], parameters[[:kf]], (i, p) -> p[:kf] * i[:F2]),
        SimpleFlux([:F3], [:Qf3], parameters[[:kf]], (i, p) -> p[:kf] * i[:F3]),
        SimpleFlux([:S1], [:Qf3], parameters[[:ks]], (i, p) -> p[:ks] * i[:S1]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            F1=input[:Pf] - input[:Qf1],
            F2=input[:Qf1] - input[:Qf2],
            F3=input[:Qf2] - input[:Qf3],
            S1=input[:Ps] - input[:S1]
        )
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end


function RoutingStore_XAJ(; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:InterRunoff, :BaseRunoff]),
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            InterRoutingStore=@.(input[:InterRunoff] - parameters[:ci] * input[:InterRoutingStore]),
            BaseRoutingStore=@.(input[:BaseRunoff] - parameters[:cg] * input[:BaseRoutingStore]),
        )
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end