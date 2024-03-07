"""
SoilWaterReservoir in Exp-Hydro
"""
function Soil_ExpHydro(; name::Symbol, parameters::ComponentVector{T}, init_states::ComponentVector{T}) where {T<:Number}
    funcs = [
        Tranparent([:Infiltration]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Baseflow([:SoilWater], parameters=parameters[[:Smax, :Qmax, :f]]),
        Surfaceflow([:SoilWater], parameters=parameters[[:Smax]]),
        Flow([:Baseflow, :Surfaceflow])
    ]

    d_funcs = [
        D_Soilwater(Dict(:In => :Infiltration, :Out => [:Evap, :Flow]), [:SoilWater])
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
function Soil_M50(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

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
        D_Soilwater([:SoilWater, :Infiltration, :Lday, :Evap, :Flow])
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
SoilWaterReservoir in M100
"""
function Soil_M100(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    ann = Lux.Chain(Lux.Dense(3, 16, tanh), Lux.Dense(16, 16, leakyrelu), Lux.Dense(16, 1, leakyrelu))

    funcs = [
        Tranparent([:Melt, :Rainfall, :Temp, :SnowWater, :SoilWater, :Evap, :Lday]),
        NNFlux([:SnowWater, :SoilWater, :Temp], [:Evap], model=et_ann, seed=42)
        #TODO 还没写完
    ]

    d_funcs = [
        D_Soilwater([:SnowWater, :SoilWater, :Rainfall, :Melt, :Lday, :Evap, :Flow])
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
SoilWaterReservoir in GR4J
"""
function Soil_GR4J(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:SoilWater, :Infiltration], parameters=parameters[[:x1]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:x1]]),
        Percolation([:SoilWater], parameters=parameters[[:x1]]),
        HydroFlux([:Infiltration, :Percolation, :Saturation], [:TempFlow],
            ComponentVector{T}(),
            (i, p, sf) -> @.(i[:Infiltration] - i[:Saturation] + i[:Percolation])),
        Splitter([:TempFlow], [:SlowFlow, :FastFlow], parameters=ComponentVector{T}(SlowFlow=0.9, FastFlow=0.1))
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Percolation]), [:SoilWater])
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
SoilWaterReservoir in HYMOD
"""
function Soil_HyMOD(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:SoilWater, :Infiltration], parameters=parameters[[:Smax, :b]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:Smax]]),
        Splitter([:Saturation], [:FastFlow, :SlowFlow],
            parameters=ComponentVector(FastFlow=parameters[:a], SlowFlow=(T(1) .- parameters[:a])))
    ]

    d_funcs = [
        Differ(Dict(:Infiltration => :In, :Evap => :Evap, :Saturation => :Out), [:SoilWater])
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
SoilWaterReservoir in XAJ
"""
function Soil_XAJ(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        Saturation([:TensionWater, :Infiltration], parameters=parameters[[:Aim, :Wmax, :a, :b]]),
        Evap([:SoilWater, :Pet], parameters=parameters[[:c, :LM]]),
    ]

    d_funcs = [
        Differ(Dict(:In => [:Infiltration], :Out => [:Evap, :Saturation]), [:SoilWater]),
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
HBV
"""
function Soil_HBV(; name::Symbol,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T}) where {T<:Number}

    funcs = [
        HydroFlux([:SoilWater], [:Capillary], parameters=parameters[[:cflux, :fc]],
            func=(i, p, sf) -> @.(p[:cflux] * (1 - i[:SoilWater] / p[:fc]))),
        Evap([:SoilWater, :Pet], parameters=parameters[[:lp, :fc]]),
        Recharge([:SoilWater, :Infiltration], parameters=parameters[[:fc, :β]]),
        HydroFlux([:UpperZone], [:InterFlow], parameters=parameters[[:k0, :α]],
            func=(i, p, sf) -> @.(p[:k0] * i[:UpperZone]^(1 + p[:α]))),
        HydroFlux([:LowerZone], [:BaseFlow], parameters=parameters[[:k1]],
            func=(i, p, sf) -> @.(p[:k1] * i[:LowerZone])),
        Summation([:InterFlow, :BaseFlow], [:Flow])
    ]


    d_funcs = [
        Differ(Dict(:In => [:Infiltration, :Capillary], :Out => [:Evap, :Recharge]), [:SoilWater]),
        HydroFlux([:Recharge, :Capillary, :InterFlow], [:UpperZone], parameters=parameters[[:c]],
            func=(i, p, sf) -> @.(i[:Recharge] - i[:Capillary] - i[:InterFlow] - p[:c])),
        HydroFlux([:BaseFlow], [:LowerZone], parameters=parameters[[:c]],
            func=(i, p, sf) -> @.(p[:c] - i[:BaseFlow]))
    ]

    ODEElement(
        name=name,
        parameters=parameters,
        init_states=init_states,
        funcs=funcs,
        d_funcs=d_funcs
    )
end