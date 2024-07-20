using ComponentArrays

# TODO 有个方法获取能够应对当前情况，我们可以先试试unit的参数合并，暂时不考虑其他的，我们可能可以使用unit的参数信息和init信息来解决这个合并问题
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

a1 = ComponentVector(a=1, b=(c=3, d=(e=3, f=1, d=3)), c=3)
a2 = ComponentVector(a=2, b=(c=2, d=(e=4, f=2)))

collect(merge_ca(a1, a2))
