# DeepFlex.jl 设计思路

## 这个包的目标是什么？

- 适用于概念性水文模型搭建，以及物理耦合的深度学习模型搭建
- 继承superflex和MARRMoT各自的优缺点，构建新的水文模型构建框架，具体来说：
  - 沿用superflex的构建框架，即element，node，network等
  - 沿用MARRMoT的flux function构建思路，使不同计算模块实现分类管理
  - **保证各模块的随机组合的可行性，当前计算模块可粗略分为SnowWater, SoilWater和RoutingStore三种**
- 提供区域多模型按权重预测（Node），提供半分布式水文模型搭建（Network）

## 这个包的使用风格应该是什么样子的？

`flux -> Element -> Unit -> Node -> Network`

- Function根据参数不同分派不同函数计算

  ```
  baseflow(i::NamedTuple, p::NamedTuple, sf::Function)
  ```
- 让概念性水文模型的搭建像torch，lux这样，每个中间通量的计算都应该是一个层

  ```julia
  model = build_unit(
      name=name,
      elements=elements
  )
  ```
- 每个element又是由多个flux function所组成的，包中已经内置了大量的flux function和构成的element

  ```
  func1 = flux_func(input...;params...)
  element = Element(
  		func1,
  		func2;
  		params...,
  		config...
  )
  ```
- 构建模型中，应该有构建合理性的校验和提示功能，主要的提示功能就是对模型各层计算数据是否存在缺失的功能，并对于缺失的模块提出增加建议
- 深度学习模型嵌入应该是作为特殊的function进行嵌入
- Element是由多个function组合而成，要保证能使用既有的func也能使用自定义的func
- 一般而言模型自身不携带任何数据,一般通过外界参数输入或参数估计器输入

## SoilWater的一些共性
- 水文模型的核心计算模块，通常代表模型的基本输入
- 模型输出为flow，包括：
  1. 单一flow，如Exphydro
  2. Fastflow和Slowflow, 如GR4J
  3. Surfaceflow, Interflow, Baseflow三层
- 模型输入需要统一，比如Infiltration，Pet