"""
RoutingStore in GR4J
"""
function Routing_ExpHydro(; name::Symbol)

    funcs = [
        FlowFlux([:baseflow, :surfaceflow])
    ]

    HydroElement(
        name=name,
        funcs=funcs
    )
end


"""
RoutingStore in GR4J
"""
function Routing_GR4J(; name::Symbol)

    funcs = [
        RechargeFlux([:routingstore], param_names=[:x2, :x3, :ω]),
        SimpleFlux([:routingstore], :routedflow,
            param_names=[:x3, :γ],
            func=(i, p, sf) -> @.((p[:x3]^(1 - p[:γ])) / (p[:γ] - 1) * (i[:routingstore]^p[:γ]))),
        SimpleFlux([:routedflow, :recharge, :fastflow], :flow,
            param_names=[:k1],
            func=(i, p, sf) -> @.(i[:routedflow] + i[:recharge] + i[:fastflow]))
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:slowflow, :recharge], :Out => [:routedflow]), :routingstore),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function Routing_HBV(; name::Symbol)

    funcs = [
        SimpleFlux([:upperzone], :interflow,
            param_names=[:k0, :α],
            func=(i, p, sf) -> @.(p[:k0] * i[:upperzone]^(1 + p[:α]))),
        SimpleFlux([:lowerzone], :baseflow,
            param_names=[:k1],
            func=(i, p, sf) -> @.(p[:k1] * i[:lowerzone])),
        SimpleFlux([:interflow, :baseflow], :flow,
            param_names=[:k1],
            func=(i, p, sf) -> @.(i[:interflow] + i[:baseflow])),
    ]

    dfuncs = [
        SimpleFlux([:recharge, :capillary, :interflow], :upperzone,
            param_names=[:c],
            func=(i, p, sf) -> @.(i[:recharge] - i[:capillary] - i[:interflow] - p[:c])),
        SimpleFlux([:baseflow], :lowerzone,
            param_names=[:c],
            func=(i, p, sf) -> @.(p[:c] - i[:baseflow]))
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

function Routing_HyMOD(; name::Symbol)

    funcs = [
        SimpleFlux([:fastrouting1], :qf1, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:fastrouting1]),
        SimpleFlux([:fastrouting2], :qf2, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:fastrouting2]),
        SimpleFlux([:fastrouting3], :qf3, param_names=[:kf], func=(i, p) -> p[:kf] .* i[:fastrouting3]),
        SimpleFlux([:slowrouting], :Qs, param_names=[:ks], func=(i, p) -> p[:ks] .* i[:slowrouting]),
        Summation([:Qs, :qf3], :flow)
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:fastflow], :Out => [:qf1]), :fastrouting1),
        DifferFlux(Dict(:In => [:qf1], :Out => [:qf2]), :fastrouting2),
        DifferFlux(Dict(:In => [:qf2], :Out => [:qf3]), :fastrouting3),
        DifferFlux(Dict(:In => [:slowflow], :Out => [:Qs]), :slowrouting),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


function Routing_XAJ(; name::Symbol)

    tmp_func = (i, p, sf) -> begin
        free_water, flux_in = i[:freewater], i[:fluxin]
        Smax, ex = p[:Smax], p[:ex]
        tmp_re = @.(sf(1 - free_water / Smax) * (1 - free_water / Smax))
        @.(1 - (sf(1 - tmp_re) * (1 - tmp_re) + sf(tmp_re - 1) * (tmp_re - 1))^ex * flux_in)
    end

    funcs = [
        SimpleFlux(Dict(:freewater => :freewater, :saturation => :fluxin), :surfacerunoff, param_names=[:Smax, :ex], func=tmp_func),
        SimpleFlux(Dict(:freewater => :freewater, :surfaceflow => :fluxin), :interrunoff, param_names=[:Smax, :ex], func=tmp_func),
        SimpleFlux(Dict(:freewater => :freewater, :interflow => :fluxin), :baserunoff, param_names=[:Smax, :ex], func=tmp_func), SimpleFlux([:interrouting], :interflow, param_names=[:ci], func=(i, p, sf) -> p[:ci] .* i[:interrouting]),
        SimpleFlux([:baserouting], :baseflow, param_names=[:cg], func=(i, p, sf) -> p[:cg] .* i[:baserouting]),
        SimpleFlux([:prcp, :surfacerunoff], :surfaceflow, param_names=[:Aim], func=(i, p, sf) -> p[:Aim] .* i[:prcp] .+ i[:surfaceflow]),
        Summation([:baseflow, :interflow, :surfaceflow], :flow)
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:saturation], :Out => [:baseflow, :interflow, :surfaceflow]), [:freewater]),
        DifferFlux(Dict(:In => [:interrunoff], :Out => [:interflow]), [:interrouting]),
        DifferFlux(Dict(:In => [:baserunoff], :Out => [:baseflow]), [:baserouting]),
    ]

    HydroElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end