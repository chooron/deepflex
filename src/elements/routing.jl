"""
RoutingStore in GR4J
"""
function Routing_ExpHydro(; name::Symbol)

    funcs = [
        Flow([:BaseFlow, :SurfaceFlow])
    ]

    SimpleElement(
        name=name,
        funcs=funcs
    )
end


"""
RoutingStore in GR4J
"""
function Routing_GR4J(; name::Symbol)

    funcs = [
        Tranparent([:SlowFlow]),
        Recharge([:RoutingStore], param_names=[:x2, :x3, :ω]),
        SimpleFlux([:RoutingStore], :RoutedFlow,
            param_names=[:x3, :γ],
            func=(i, p, sf) -> @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:RoutingStore]^p[:γ]))),
        Summation([:RoutedFlow, :Recharge, :FastFlow], :Flow)
    ]

    d_funcs = [
        Differ(Dict(:In => [:SlowFlow, :Recharge], :Out => [:RoutedFlow]), [:RoutingStore]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

function Routing_HBV(; name::Symbol)

    funcs = [
        Tranparent([:Recharge, :Capillary]),
        SimpleFlux([:UpperZone], :InterFlow,
            param_names=[:k0, :α],
            func=(i, p, sf) -> @.(p[:k0] * i[:UpperZone]^(1 + p[:α]))),
        SimpleFlux([:LowerZone], :BaseFlow,
            param_names=[:k1],
            func=(i, p, sf) -> @.(p[:k1] * i[:LowerZone])),
        Summation([:InterFlow, :BaseFlow], :Flow)
    ]

    d_funcs = [
        SimpleFlux([:Recharge, :Capillary, :InterFlow], :UpperZone,
            param_names=[:c],
            func=(i, p, sf) -> @.(i[:Recharge] - i[:Capillary] - i[:InterFlow] - p[:c])),
        SimpleFlux([:BaseFlow], :LowerZone,
            param_names=[:c],
            func=(i, p, sf) -> @.(p[:c] - i[:BaseFlow]))
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

function Routing_HyMOD(; name::Symbol)

    funcs = [
        Tranparent([:FastFlow, :SlowFlow]),
        SimpleFlux([:FastRouting1], :Qf1, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:FastRouting1]),
        SimpleFlux([:FastRouting2], :Qf2, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:FastRouting2]),
        SimpleFlux([:FastRouting3], :Qf3, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:FastRouting3]),
        SimpleFlux([:SlowRouting], :Qs, param_names=[:ks], func=(i, p) -> p[:ks] .* i[:SlowRouting]),
        Summation([:Qs, :Qf3], :Flow)
    ]

    d_funcs = [
        Differ(Dict(:In => [:FastFlow], :Out => [:Qf1]), [:FastRouting1]),
        Differ(Dict(:In => [:Qf1], :Out => [:Qf2]), [:FastRouting2]),
        Differ(Dict(:In => [:Qf2], :Out => [:Qf3]), [:FastRouting3]),
        Differ(Dict(:In => [:SlowFlow], :Out => [:Qs]), [:SlowRouting]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


function Routing_XAJ(; name::Symbol)

    tmp_func = (i, p, sf) -> begin
        free_water, flux_in = i[:FreeWater], i[:FluxIn]
        Smax, ex = p[:Smax], p[:ex]
        tmp_re = @.(sf(1 - free_water / Smax) * (1 - free_water / Smax))
        @.(1 - (sf(1 - tmp_re) * (1 - tmp_re) + sf(tmp_re - 1) * (tmp_re - 1))^ex * flux_in)
    end

    funcs = [
        SimpleFlux(Dict(:FreeWater => :FreeWater, :Saturation => :FluxIn), :SurfaceRunoff, param_names=[:Smax, :ex], func=tmp_func),
        SimpleFlux(Dict(:FreeWater => :FreeWater, :SurfaceFlow => :FluxIn), :InterRunoff, param_names=[:Smax, :ex], func=tmp_func),
        SimpleFlux(Dict(:FreeWater => :FreeWater, :InterFlow => :FluxIn), :BaseRunoff, param_names=[:Smax, :ex], func=tmp_func),
        
        SimpleFlux([:InterRouting], :InterFlow, param_names=[:ci], func=(i, p, sf) -> p[:ci] .* i[:InterRouting]),
        SimpleFlux([:BaseRouting], :BaseFlow, param_names=[:cg], func=(i, p, sf) -> p[:cg] .* i[:BaseRouting]),
        SimpleFlux([:Prcp, :SurfaceRunoff], :SurfaceFlow, param_names=[:Aim], func=(i, p, sf) -> p[:Aim] .* i[:Prcp] .+ i[:SurfaceFlow]),
        Summation([:BaseFlow, :InterFlow, :SurfaceFlow], :Flow)
    ]

    d_funcs = [
        Differ(Dict(:In => [:Saturation], :Out => [:BaseFlow, :InterFlow, :SurfaceFlow]), [:FreeWater]),
        Differ(Dict(:In => [:InterRunoff], :Out => [:InterFlow]), [:InterRouting]),
        Differ(Dict(:In => [:BaseRunoff], :Out => [:BaseFlow]), [:BaseRouting]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end