using ComponentArrays

ca1 = ComponentVector(a=1, b=[2, 3, 4], c=(b=1, c=2))
ax1 = getaxes(ca1)

ComponentVector(ones(length(ca1)),ax1)