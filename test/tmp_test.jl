using ComponentArrays

# 假设有一个 ComponentVector 的结果
re = ComponentVector(
    SnowWater=0.0, SoilWater=1303.004248)
result = [re, re, re]
# 将相同名称的组件拼接成一个 Vector
merged_vector = hcat(result...)

# 打印结果
println(merged_vector)
