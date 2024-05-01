function merge_ca(ca1::ComponentArray, ca2::ComponentArray, key::Symbol=:param)
    #* used for merge data
    share_keys = intersect(keys(ca1), keys(ca2))
    # println((eltype(ca1)))
    # println((eltype(ca2)))
    new_ca1 = ComponentVector{promote_type(eltype(ca1),eltype(ca2))}() # {eltype(ca1)}
    if length(share_keys) > 0
        for key in share_keys
            if !(typeof(ca1[key]) <: ComponentArray) | !(typeof(ca2[key]) <: ComponentArray)
                new_ca1 = ComponentVector(new_ca1; ComponentVector([key=>ca2[key]])...)
            else
                merged_ca = merge_ca(ca1[key], ca2[key], key)
                new_ca1 = ComponentVector(new_ca1; merged_ca...)
            end
        end
    end
    for key in filter(k -> !(k in share_keys), keys(ca1))
        new_ca1 = ComponentVector(new_ca1; ca1[[key]]...)
    end
    for key in filter(k -> !(k in share_keys), keys(ca2))
        new_ca1 = ComponentVector(new_ca1; ca2[[key]]...)
    end
    ComponentVector(namedtuple([key],[new_ca1]))
end