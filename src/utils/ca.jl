
# merge ComponentArray:https://github.com/jonniedie/ComponentArrays.jl/issues/186
function merge_ca(ca::ComponentArray{T1}, ca2::ComponentArray{T2}) where {T1,T2}
    ax = getaxes(ca)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    _p = Vector{T2}()
    for vk in vks
        if length(getaxes(ca[vk])) > 0
            _p = vcat(_p, collect(merge_ca(ca[vk], vk in vks2 ? getproperty(ca2, vk) : ComponentVector())))
        else
            if vk in vks2
                _p = vcat(_p, ca2[vk])
            else
                _p = vcat(_p, ca[vk])
            end
        end
    end
    ComponentArray(_p, ax)
end