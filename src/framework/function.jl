
function get_parameters(func::Func; names::Vector{Symbol}=nothing)::Dict{Symbol,Any}
    if isnothing(names)
        return func.parameters
    else
        return Dict(name => func.parameters[name] for name in names)
    end
end

function set_parameters!(func::Func; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        if p.name in keys(func.parameters)
            func.parameters[p.name] = p.value
        end
    end
end