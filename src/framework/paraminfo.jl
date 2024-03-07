mutable struct ParamInfo{T} <: AbstractParamInfo where {T<:Number}
    name::Symbol
    default::T
    lb::T
    ub::T
    value::T
end

function ParamInfo(name::Symbol, default::T; lb::T=nothing, ub::T=nothing) where {T<:Number}
    lb = isnothing(lb) ? default : lb
    ub = isnothing(ub) ? default : ub
    ParamInfo{T}(name, default, lb, ub, default)
end

function update_paraminfo!(paraminfo::P, paramvalue::T) where {P<:AbstractParamInfo,T<:Number}
    paraminfo.value = paramvalue
end

function update_paraminfos!(paraminfos::Vector{P}, paramvalues::Vector{T}) where {P<:AbstractParamInfo,T<:Number}
    for (idx, value) in enumerate(paramvalues)
        update_paraminfo!(paraminfos[idx], value)
    end
end

const NO_PARAMETER = NamedTuple()