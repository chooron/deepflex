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
        HydroFlux([:RoutingStore], [:RoutedFlow], parameters=parameters[[:x3, :γ]],
            func=(i, p, sf) -> @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:RoutingStore]^p[:γ]))),
        Summation([:RoutedFlow, :Recharge, :FastFlow], [:Flow])
    ]

    d_funcs = [
        Differ(Dict(:In => [:SlowFlow, :Recharge], :Out => [:RoutedFlow]), [:RoutingStore]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

function Routing_HyMOD(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Tranparent([:FastFlow, :SlowFlow]),
        SimpleFlux([:FastRouting1], [:Qf1], parameters[[:kf]], (i, p) -> [p[:kf] .* i[:FastRouting1]]),
        SimpleFlux([:FastRouting2], [:Qf2], parameters[[:kf]], (i, p) -> [p[:kf] .* i[:FastRouting2]]),
        SimpleFlux([:FastRouting3], [:Qf3], parameters[[:kf]], (i, p) -> [p[:kf] .* i[:FastRouting3]]),
        SimpleFlux([:SlowRouting], [:Qs], parameters[[:ks]], (i, p) -> [p[:ks] .* i[:SlowRouting]]),
        Summation([:Qs, :Qf3], [:Flow])
    ]

    d_funcs = [
        Differ(Dict(:In => [:FastFlow], :Out => [:Qf1]), [:FastRouting1]),
        Differ(Dict(:In => [:Qf1], :Out => [:Qf2]), [:FastRouting2]),
        Differ(Dict(:In => [:Qf2], :Out => [:Qf3]), [:FastRouting3]),
        Differ(Dict(:In => [:SlowFlow], :Out => [:Qs]), [:SlowRouting]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


function Routing_XAJ(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    tmp_func = (i, p, sf) -> begin
        free_water, flux_in = i[:FreeWater], i[:FluxIn]
        Smax, ex = p[:Smax], p[:ex]
        tmp_re = @.(sf(1 - free_water / Smax) * (1 - free_water / Smax))
        @.(1 - (sf(1 - tmp_re) * (1 - tmp_re) + sf(tmp_re - 1) * (tmp_re - 1))^ex * flux_in)
    end

    funcs = [
        HydroFlux(Dict(:FreeWater => :FreeWater, :Saturation => :FluxIn), [:SurfaceRunoff],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
        HydroFlux(Dict(:FreeWater => :FreeWater, :SurfaceFlow => :FluxIn), [:InterRunoff],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
        HydroFlux(Dict(:FreeWater => :FreeWater, :InterFlow => :FluxIn), [:BaseRunoff],
            parameters=parameters[[:Smax, :ex]], func=tmp_func),
        HydroFlux([:InterRouting], [:InterFlow],
            parameters=parameters[:ci],
            func=(i, p, sf) -> p[:ci] .* i[:InterRouting]),
        HydroFlux([:BaseRouting], [:BaseFlow],
            parameters=parameters[:cg],
            func=(i, p, sf) -> p[:cg] .* i[:BaseRouting]),
        HydroFlux([:Prcp, :SurfaceRunoff], [:SurfaceFlow],
            parameters=parameters[[:Aim]],
            func=(i, p, sf) -> p[:Aim] .* i[:Prcp] .+ i[:SurfaceFlow]),
        Summation([:BaseFlow, :InterFlow, :SurfaceFlow], [:Flow])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Saturation], :Out => [:BaseFlow, :InterFlow, :SurfaceFlow]), [:FreeWater]),
        Differ(Dict(:In => [:InterRunoff], :Out => [:InterFlow]), [:InterRouting]),
        Differ(Dict(:In => [:BaseRunoff], :Out => [:BaseFlow]), [:BaseRouting]),
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        d_funcs=d_funcs
    )
end