"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil_ExpHydro(; name::Symbol)
    funcs = [
        Tranparent([:Infiltration]),
        Evap([:SoilWater, :Pet], parameter_names=[:Smax]),
        Baseflow([:SoilWater], parameter_names=[:Smax, :Qmax, :f]),
        Surfaceflow([:SoilWater], parameter_names=[:Smax]),
        Flow([:Baseflow, :Surfaceflow])
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Flow]), [:SoilWater])
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
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
            parameter_names=Symbol[],
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
        Tranparent([:Melt, :Rainfall, :Temp, :SnowWater, :SoilWater, :Evap, :Lday]),
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=ann, seed=42)
        #TODO 还没写完
    ]

    d_funcs = [
        D_Soilwater([:SnowWater, :SoilWater, :Rainfall, :Melt, :Lday, :Evap, :Flow])
    ]

    d_funcs = [
        SimpleFlux([:SnowWater, :SoilWater, :Rainfall, :Melt, :Lday, :Evap, :Flow], [:SoilWater],
            parameter_names=Symbol[],
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
        Saturation([:SoilWater, :Infiltration], parameter_names=[:x1]),
        Evap([:SoilWater, :Pet], parameter_names=[:x1]),
        Percolation([:SoilWater], parameter_names=[:x1]),
        SimpleFlux([:Infiltration, :Percolation, :Saturation], :TempFlow,
            parameters=Symbol[],
            func=(i, p, sf) -> @.(i[:Infiltration] - i[:Saturation] + i[:Percolation])),
        Splitter([:TempFlow], [:SlowFlow, :FastFlow], parameters=ComponentVector{T}(SlowFlow=0.9, FastFlow=0.1))
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
        Saturation([:SoilWater, :Infiltration], parameter_names=[:Smax, :b]),
        Evap([:SoilWater, :Pet], parameter_names=[:Smax]),
        SimpleFlux([:Saturation], [:FastFlow, :SlowFlow], parameter_names=[:a],
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
        Saturation([:SoilWater, :Infiltration], parameter_names=[:Aim, :Wmax, :a, :b]),
        Evap([:SoilWater, :Pet], parameter_names=[:c, :LM]),
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
            parameter_names=[:cflux, :fc],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:SoilWater] / p[:fc]))),
        Evap([:SoilWater, :Pet], parameter_names=[:lp, :fc]),
        Recharge([:SoilWater, :Infiltration], parameter_names=[:fc, :β]),
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