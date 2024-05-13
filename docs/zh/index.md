```@meta
CurrentModule = LumpedHydro
```

# LumpedHydro.jl

LumpedHydro.jl是一个基于julia语言编写的用于构建概念性水文模型的包，通过这个包可以构建从单一的计算模块到一个完整的概念性水文模型再到一个河网式分布式水文模型。

## 特性

* [X] 基于julia语言面向函数、多重分派的编程特性，形成了名称驱动的模型构建方式（见自定义模型构建）

  * [X] 通过输入输出以及参数名称快速搜索某一个通量的计算函数（见名称驱动的函数搜索方式）
  * [X] 获取模块或模型的输入输出以及参数名称自动生成模型运行/优化所需参数的规范
  * [X] 但是在一些情况中可能会创建一些并没有实际含义的变量
* [X] 构建no MTK和MTK风格(MTK风格参见ModelingToolkit.jl)的常微分方程求解方法，适应不同需求的模型应用需求

  * [X] no MTK风格通常确保最基础的计算能力，能够保证复杂情况的运行，但是计算性能略低于MTK风格
  * [X] MTK风格基于ModelingToolkit.jl构建ODESystem，有着更高的计算性能，但未知错误较多（见性能比较）
* [X] 统一式与逐步式的多模块或模型求解方式

  * [X] 统一式求解方式是将多个模块中的通量计算函数整合到一个常微分方程中，减小中间变量的插值损失与计算成本
  * [X] 逐步式求解方式则是针对不同模块进行独立求解
* [X] 基于网络拓扑计算的通量函数自动匹配计算和河网汇流的网络逻辑计算

  * [X] 通量函数之间可以通过输入输出的变量名称构建计算网，在搭建模块时能够不考虑通量函数的顺序，支撑模块临时添加通量函数的功能
  * [X] 河网汇流计算是基于各个节点的拓扑连接关系实现，需要并行计算各个节点的内置计算过程后按河网网络顺序计算实现各子流域汇流计算
* [X] 提供NeuralFlux对象支持神经网络模型耦合水文模型，主要包括计算公式替换，模型参数估计

  * [X] 在计算公式替换中，这个包可以根据Lux.jl构建的神经网络构建NeuralFlux并用于计算水文模型中某些通量，典型模型由m50
  * [X] 在模型参数估计中，这个包将参数作为一种无量纲的通量，并根据NeuralFlux作出预测，典型的有dPL模型
* [X] 基于Optimization.jl提供模型参数的全局搜索（黑箱）和局部搜索（梯度）的优化算法

  * [X] 提供BoxOptimization或更多算法，用于初步的模型参数优化
  * [X] 局部搜索或梯度则更多的用于神经网络内部权重的优化

## 安装

To install ModelingToolkitNeuralNets.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("LumpedHydro")
```

## 贡献

- Please refer to the
  [SciML ColPrac: Contributor&#39;s Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to SciML.
- See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
- There are a few community forums:

  + The #diffeq-bridged and #sciml-bridged channels in the
    [Julia Slack](https://julialang.org/slack/)
  + The #diffeq-bridged and #sciml-bridged channels in the
    [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
  + On the [Julia Discourse forums](https://discourse.julialang.org)
  + See also [SciML Community page](https://sciml.ai/community/)
