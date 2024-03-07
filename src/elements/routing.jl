"""
RoutingStore in GR4J
"""
function Routing_ExpHydro(; name::Symbol)

    funcs = [
        Flow([:Baseflow, :Surfaceflow])
    ]

    SimpleElement(
        name=name,
        parameters=ComponentVector(),
        funcs=funcs
    )
end


"""
RoutingStore in GR4J
"""
function Routing_GR4J(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:SlowFlow]),
        Recharge([:RoutingStore], parameters=parameters[[:x2, :x3, :ω]]),
        SimpleFlux([:RoutingStore], [:RoutedFlow], parameters=parameters[[:x3, :γ]],
            func=(i, p) -> [@.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:RoutingStore]^p[:γ]))]),
        Summation([:RoutedFlow, :Recharge, :FastFlow], [:Flow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(RoutingStore=input[:SlowFlow] + input[:Recharge] - input[:RoutedFlow])
    end

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        get_du=get_du
    )
end

function Routing_HyMOD(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:FastFlow, :SlowFlow]),
        SimpleFlux([:FastRouting1], [:Qf1], parameters[[:kf]], (i, p) -> [p[:kf] * i[:FastRouting1]]),
        SimpleFlux([:FastRouting2], [:Qf2], parameters[[:kf]], (i, p) -> [p[:kf] * i[:FastRouting2]]),
        SimpleFlux([:FastRouting3], [:Qf3], parameters[[:kf]], (i, p) -> [p[:kf] * i[:FastRouting3]]),
        SimpleFlux([:SlowRouting], [:Qs], parameters[[:ks]], (i, p) -> [p[:ks] * i[:SlowRouting]]),
        Summation([:Qs, :Qf3], [:Flow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            FastRouting1=input[:FastFlow] - input[:Qf1],
            FastRouting2=input[:Qf1] - input[:Qf2],
            FastRouting3=input[:Qf2] - input[:Qf3],
            SlowRouting=input[:SlowFlow] - input[:Qs]
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


function Routing_XAJ(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:BaseFlow, :InterFlow, :SurfaceFlow], [:BaseRunoff, :InterRunoff, :SurfaceRunoff]),
        SimpleFlux([:InterRouting], [:InterFlow], parameters=ComponentVector(p=parameters[:ci]), func=(i, p) -> [p[:p] .* i[:InterRouting]]),
        SimpleFlux([:BaseRouting], [:BaseFlow], parameters=ComponentVector(p=parameters[:cg]), func=(i, p) -> [p[:p] .* i[:BaseRouting]]),
        SimpleFlux([:Prcp, :SurfaceRunoff], [:SurfaceFlow], parameters=parameters[[:Aim]], func=(i, p) -> [p[:Aim] .* i[:Prcp] .+ i[:SurfaceFlow]]),
        Summation([:BaseFlow, :InterFlow, :SurfaceFlow], [:Flow])
    ]

    get_du = (input::ComponentVector{T}, parameters::ComponentVector{T}) -> begin
        ComponentVector(
            InterRouting=@.(input[:InterRunoff] - input[:InterFlow]),
            BaseRouting=@.(input[:BaseRunoff] - input[:BaseFlow]),
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