using ComponentArrays

# TODO 有个方法获取能够应对当前情况，我们可以先试试unit的参数合并，暂时不考虑其他的，我们可能可以使用unit的参数信息和init信息来解决这个合并问题
function merge_ca(ca::ComponentArray{T1}, ca2::ComponentArray{T2}) where {T1,T2}
    ax = getaxes(ca)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    _p = Vector{T2}()
    for vk in vks
        if vk in vks2
            _p = vcat(_p, ca2[vk])
        else
            _p = vcat(_p, ca[vk])
        end
    end
    ComponentArray(_p, ax)
end

a1 = ComponentVector(a=1, b=(c=3, d=4), c=3)
a2 = ComponentVector(a=2, b=(c=2, d=5))

merge_ca(a1, a2)
 