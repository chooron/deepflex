```@meta
CurrentModule = LumpedHydro
```

# LumpedHydro.jl

LumpedHydro.jl是一个基于julia语言编写的用于构建概念性水文模型的包，通过这个包可以构建从单一的计算模块到一个完整的概念性水文模型。

## 安装

To install ModelingToolkitNeuralNets.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("LumpedHydro")
```

## 运行一个ExpHydro水文模型

```julia
# import lib
using CSV
using DataFrames
using CairoMakie
using BenchmarkTools
using ComponentArrays
using LumpedHydro

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
flow_vec = df[ts, "flow(mm)"]

# build model
model = LumpedHydro.ExpHydro.Unit(name=:exphydro, mtk=false)

# setup parameters and initstates
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
unit_params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_init_states = (snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector((params=unit_params, initstates=unit_init_states, weight=1.0))

# run model
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
solver = LumpedHydro.ODESolver(reltol=1e-3, abstol=1e-3)
result = model(input, pas, timeidx=ts, solver=solver);
result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts, flow_vec, color=:red)
lines!(ax, ts, result_df[!,"flow"], color=:blue)
fig
```

![predictions](docs/picture/predictions.png)

## 特性

* [X] 基于julia语言面向函数、多重分派的编程特性，形成了名称驱动的模型构建方式（见自定义模型构建）

  * [X] 通过输入输出以及参数名称快速搜索某一个通量的计算函数（见名称驱动的函数搜索方式）
  * [X] 获取模块或模型的输入输出以及参数名称自动生成模型运行/优化所需参数的规范
  * [X] 能够实现在无实际参数情况下完成模型的构建，实现模型与模型参数的完全解耦
  * [X] 但是在一些情况中可能会创建一些并没有实际含义的变量

* [X] 提供NeuralFlux对象支持神经网络模型耦合水文模型，主要包括计算公式替换，模型参数估计

  * [X] 在计算公式替换中，这个包可以根据Lux.jl构建的神经网络构建NeuralFlux并用于计算水文模型中某些通量，典型模型由m50
  * [X] 在模型参数估计中，这个包将参数作为一种无量纲的通量，并根据NeuralFlux作出预测，典型的有dPL模型
* [X] 基于Optimization.jl提供模型参数的全局搜索（黑箱）和局部搜索（梯度）的优化算法

  * [X] 提供BoxOptimization或更多算法，用于初步的模型参数优化
  * [X] 局部搜索或梯度则更多的用于神经网络内部权重的优化
* [X] 构建No MTK和MTK风格(MTK风格参见ModelingToolkit.jl)的常微分方程求解方法，适应不同需求的模型应用需求

  * [X] No MTK风格通常确保最基础的计算能力，能够保证复杂情况的运行，但是计算性能略低于MTK风格
  * [X] MTK风格基于ModelingToolkit.jl构建ODESystem，有着更高的计算性能，但未知错误较多（见性能比较）

## 代码存在的不足

- Zygote不支持MTK存在的try-catch（通常是在参数设置中引起的）: `ERROR: Compiling Tuple{Type{Dict}, Dict{Any, Any}}: try/catch is not supported.`同时Zygote也不支持mutable array，因此也不支持非MTK写法，而Enzyme需要先对代码进行编译，因此存在较多类型问题。
- MTK 与zygote的兼容性实在是太差了，包括数据输入设置，参数设定等

## 未来工作计划

* [ ] 复现当前主流的概念水文模型和已发表并具有代表性的神经耦合的水文模型
* [ ] 从Lux.jl中LSTM Layer实现方式，来借鉴并修改非MTK的写法，
* [ ] 根据构建的网络使用迭代计算方式不断从上层调度函数，以此避免计算资源的损耗
