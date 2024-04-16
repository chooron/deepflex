using ComponentArrays
using NamedTupleTools

function merge_ca(ca1::ComponentArray, ca2::ComponentArray)
    # 这里需要写一个名字相同时合并的方法
    share_keys = intersect(keys(ca1), keys(ca2))
    new_ca1 = ComponentVector{eltype(ca1)}()
    if length(share_keys) > 0
        for key in share_keys
            merged_ca = merge_by_key(ca1, ca2, key)
            new_ca1 = ComponentVector(new_ca1; merged_ca...)
        end
    else
        new_ca1 = ComponentVector(ca1; ca2...)
    end
    new_ca1
end

function merge_by_key(ca1::ComponentArray, ca2::ComponentArray, key::Symbol)
    # 这里需要写一个名字相同时合并的方法
    share_keys = intersect(keys(ca1[key]), keys(ca2[key]))
    if length(share_keys) > 0
        for sk in share_keys
            merge_by_key(ca1[key], ca2[key], sk)
        end
    else
        ca1 = ComponentVector(namedtuple([key], [ComponentVector(ca1[key]; ca2[key]...)]))
    end
    return ca1
end

data1 = ComponentVector(u1=(a=1, b=2))
data2 = ComponentVector(u1=(c=1, d=2))

new_data1 = merge_ca(data1, data2)