"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil_ExpHydro(; name::Symbol)
    funcs = [
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        BaseflowFlux([:soilwater], param_names=[:Smax, :Qmax, :f]),
        SurfaceflowFlux([:soilwater], param_names=[:Smax]),
        FlowFlux([:baseflow, :surfaceflow])
    ]

    dfuncs = [
        Differ(Dict(:In => [:infiltration], :Out => [:evap, :flow]), :soilwater)
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


"""
SoilWaterReservoir in M50
"""
function Soil_M50(; name::Symbol)

    et_ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
    q_ann = Lux.Chain(Lux.Dense(2, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        # ET ANN
        NNFlux([:SnowWater, :soilwater, :Temp], [:evap], model=et_ann, seed=42),
        # Q ANN
        NNFlux([:soilwater, :Prcp], [:flow], model=q_ann, seed=42),
    ]

    dfuncs = [
        SimpleFlux([:soilwater, :infiltration, :lday, :evap, :flow], :soilwater,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(input[:infiltration] -
                                  sf(input[:soilwater]) * input[:lday] * exp(input[:evap]) -
                                  sf(input[:soilwater]) * exp(input[:flow])))
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SoilWaterReservoir in M100
"""
function Soil_M100(; name::Symbol)

    ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        NNFlux([:SnowWater, :soilwater, :Temp], [:evap], model=ann, seed=42)
        #TODO 还没写完
    ]

    dfuncs = [
        SimpleFlux([:SnowWater, :soilwater, :Rainfall, :Melt, :lday, :evap, :flow], :soilwater,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(relu(sinh(input[:Rainfall])) +
                                  relu(step_func(input[:SnowWater]) * sinh(input[:Melt])) -
                                  step_func(input[:soilwater]) * input[:lday] * exp(input[:evap]) -
                                  step_func(input[:soilwater]) * exp(input[:flow])))
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SoilWaterReservoir in GR4J
"""
function Soil_GR4J(; name::Symbol)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:x1]),
        EvapFlux([:soilwater, :pet], param_names=[:x1]),
        PercolationFlux([:soilwater], param_names=[:x1]),
        SimpleFlux([:infiltration, :percolation, :saturation], :tempflow,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:infiltration] - i[:saturation] + i[:percolation])),
        SimpleFlux([:tempflow], [:slowflow, :fastflow],
            param_names=Symbol[],
            func=(i, p, sf) -> @.[i[:tempflow] * 0.9, i[:tempflow] * 0.1])
    ]

    dfuncs = [
        Differ(Dict(:In => [:infiltration], :Out => [:evap, :percolation]), :soilwater)
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function Soil_HyMOD(; name::Symbol)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:Smax, :b]),
        EvapFlux([:soilwater, :pet], param_names=[:Smax]),
        SimpleFlux([:saturation], [:fastflow, :slowflow], param_names=[:a],
            func=(i, p, sf) -> @.[i[:saturation] * (1 - p[:a]), i[:saturation] * p[:a]])
    ]

    dfuncs = [
        DifferFlux(Dict(:infiltration => :In, :evap => :evap, :saturation => :Out), :soilwater)
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


"""
SoilWaterReservoir in XAJ
"""
function Soil_XAJ(; name::Symbol)

    funcs = [
        SaturationFlux([:soilwater, :infiltration], param_names=[:Aim, :Wmax, :a, :b]),
        EvapFlux([:soilwater, :pet], param_names=[:c, :LM]),
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:infiltration], :Out => [:evap, :saturation]), :soilwater),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


"""
HBV
"""
function Soil_HBV(; name::Symbol)

    funcs = [
        SimpleFlux([:soilwater], :capillary,
            param_names=[:cflux, :fc],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:soilwater] / p[:fc]))),
        EvapFlux([:soilwater, :pet], param_names=[:lp, :fc]),
        RechargeFlux([:soilwater, :infiltration], param_names=[:fc, :β]),
    ]


    dfuncs = [
        DifferFlux(Dict(:In => [:infiltration, :capillary], :Out => [:evap, :recharge]), :soilwater),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end