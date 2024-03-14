"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil_ExpHydro(; name::Symbol)
    funcs = [
        Evap([:SoilWater, :Pet], param_names=[:Smax]),
        BaseFlow([:SoilWater], param_names=[:Smax, :Qmax, :f]),
        SurfaceFlow([:SoilWater], param_names=[:Smax]),
        Flow([:BaseFlow, :SurfaceFlow])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Flow]), [:SoilWater])
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


"""
SoilWaterReservoir in M50
"""
function Soil_M50(; name::Symbol)

    et_ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
    q_ann = Lux.Chain(Lux.Dense(2, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        Tranparent([:Infiltration]),
        # ET ANN
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=et_ann, seed=42),
        # Q ANN
        NNFlux([:SoilWater, :Prcp], [:Flow], model=q_ann, seed=42),
    ]

    d_funcs = [
        SimpleFlux([:SoilWater, :Infiltration, :Lday, :Evap, :Flow], [:SoilWater],
            param_names=Symbol[],
            func=(i, p, sf) -> @.(input[:Infiltration] -
                                  sf(input[:SoilWater]) * input[:Lday] * exp(input[:Evap]) -
                                  sf(input[:SoilWater]) * exp(input[:Flow])))
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

"""
SoilWaterReservoir in M100
"""
function Soil_M100(; name::Symbol)

    ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=ann, seed=42)
        #TODO 还没写完
    ]

    d_funcs = [
        SimpleFlux([:SnowWater, :SoilWater, :Rainfall, :Melt, :Lday, :Evap, :Flow], [:SoilWater],
            param_names=Symbol[],
            func=(i, p, sf) -> @.(relu(sinh(input[:Rainfall])) +
                                  relu(step_func(input[:SnowWater]) * sinh(input[:Melt])) -
                                  step_func(input[:SoilWater]) * input[:Lday] * exp(input[:Evap]) -
                                  step_func(input[:SoilWater]) * exp(input[:Flow])))
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

"""
SoilWaterReservoir in GR4J
"""
function Soil_GR4J(; name::Symbol)

    funcs = [
        Saturation([:SoilWater, :Infiltration], param_names=[:x1]),
        Evap([:SoilWater, :Pet], param_names=[:x1]),
        Percolation([:SoilWater], param_names=[:x1]),
        SimpleFlux([:Infiltration, :Percolation, :Saturation], :TempFlow,
            param_names=Symbol[],
            func=(i, p, sf) -> @.(i[:Infiltration] - i[:Saturation] + i[:Percolation])),
        SimpleFlux([:TempFlow], [:SlowFlow, :FastFlow],
            param_names=Symbol[],
            func=(i, p, sf) -> @.[i[:TempFlow] * 0.9, i[:TempFlow] * 0.1])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Percolation]), [:SoilWater])
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end

"""
SoilWaterReservoir in HYMOD
"""
function Soil_HyMOD(; name::Symbol)

    funcs = [
        Saturation([:SoilWater, :Infiltration], param_names=[:Smax, :b]),
        Evap([:SoilWater, :Pet], param_names=[:Smax]),
        SimpleFlux([:Saturation], [:FastFlow, :SlowFlow], param_names=[:a],
            func=(i, p, sf) -> @.[i[:Saturation] * (1 - p[:a]), i[:Saturation] * p[:a]])
    ]

    d_funcs = [
        Differ(Dict(:Infiltration => :In, :Evap => :Evap, :Saturation => :Out), [:SoilWater])
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


"""
SoilWaterReservoir in XAJ
"""
function Soil_XAJ(; name::Symbol)

    funcs = [
        Saturation([:SoilWater, :Infiltration], param_names=[:Aim, :Wmax, :a, :b]),
        Evap([:SoilWater, :Pet], param_names=[:c, :LM]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Saturation]), [:SoilWater]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end


"""
HBV
"""
function Soil_HBV(; name::Symbol)

    funcs = [
        SimpleFlux([:SoilWater], :Capillary,
            param_names=[:cflux, :fc],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:SoilWater] / p[:fc]))),
        Evap([:SoilWater, :Pet], param_names=[:lp, :fc]),
        Recharge([:SoilWater, :Infiltration], param_names=[:fc, :β]),
    ]


    d_funcs = [
        Differ(Dict(:In => [:Infiltration, :Capillary], :Out => [:Evap, :Recharge]), [:SoilWater]),
    ]

    ODEElement(
        name=name,
        funcs=funcs,
        d_funcs=d_funcs
    )
end