function SoilElement(; name::Symbol,
    funcs::Vector,
    dfuncs::Vector=SimpleFlux[],
    lfuncs::Vector=LagFlux[]
)
    # todo 针对soil element可能有着不同的判断策略

    HydroElement(
        name=Symbol(name, :_soil_),
        funcs=funcs,
        dfuncs=dfuncs,
        lfuncs=lfuncs
    )
end

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
        DifferFlux(Dict(:In => [:infiltration], :Out => [:evap, :flow]), :soilwater)
    ]

    SoilElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end


"""
SoilWaterReservoir in M50
"""
function Soil_M50(; name::Symbol)

    # 神经网络的定义是在模型之内，需要提取到模型的参数
    et_ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))
    q_ann = Lux.Chain(Lux.Dense(2, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        # normalize
        StdNormFlux(:snowwater, :norm_snowwater),
        StdNormFlux(:soilwater, :norm_soilwater),
        StdNormFlux(:temp, :norm_temp),
        StdNormFlux(:prcp, :norm_prcp),
        # ET ANN
        LuxNNFlux([:norm_snowwater, :norm_soilwater, :norm_temp], :evap, param_names=:etnn, model=et_ann, seed=42),
        # Q ANN
        LuxNNFlux([:norm_soilwater, :norm_prcp], :flow,  param_names=:qnn, model=q_ann, seed=42),
    ]

    dfuncs = [
        SimpleFlux([:soilwater, :infiltration, :lday, :evap, :flow], :soilwater, param_names=Symbol[],
            func=(i, p, sf) -> @.(input[:infiltration] -
                                  sf(input[:soilwater]) * input[:lday] * exp(input[:evap]) -
                                  sf(input[:soilwater]) * exp(input[:flow]))
        )
    ]

    SoilElement(
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

    SoilElement(
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
        SimpleFlux([:tempflow], :slowflow, param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.9),
        SimpleFlux([:tempflow], :fastflow, param_names=Symbol[], func=(i, p, sf) -> i[:tempflow] .* 0.1)
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:saturation], :Out => [:evap, :percolation]), :soilwater)
    ]

    SoilElement(
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
        SimpleFlux([:saturation], :fastflow, param_names=[:a],
            func=(i, p, sf) -> @.(i[:saturation] * (1 - p[:a]))),
        SimpleFlux([:saturation], :slowflow, param_names=[:a],
            func=(i, p, sf) -> @.(i[:saturation] * p[:a]))
    ]

    dfuncs = [
        DifferFlux(Dict(:In => [:saturation], :Out => [:evap, :saturation]), :soilwater)
    ]

    SoilElement(
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

    SoilElement(
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
        SimpleFlux([:soilwater], :capillary, param_names=[:cflux, :fc],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:soilwater] / p[:fc]))),
        EvapFlux([:soilwater, :pet], param_names=[:lp, :fc]),
        RechargeFlux([:soilwater, :infiltration], param_names=[:fc, :β]),
    ]


    dfuncs = [
        DifferFlux(Dict(:In => [:infiltration, :capillary], :Out => [:evap, :recharge]), :soilwater),
    ]

    SoilElement(
        name=name,
        funcs=funcs,
        dfuncs=dfuncs
    )
end