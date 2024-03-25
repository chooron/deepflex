using ComponentArrays

v1 = ComponentVector(a=1, b=2, c=3)
keys(v1)
@static_unpack v2 = v1;