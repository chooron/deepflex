# DeepFlex.jl 设计思路

## 这个包的目标是什么？

- 适用于概念性水文模型搭建，以及物理耦合的深度学习模型搭建
- 继承superflex和MARRMoT各自的优缺点，构建新的水文模型构建框架，具体来说：
  - 沿用superflex的构建框架，即element，node，network等
  - 沿用MARRMoT的flux function构建思路，使不同计算模块实现分类管理
  - **保证各模块的随机组合的可行性，当前计算模块可粗略分为SnowWater, SoilWater和RoutingStore三种**
- 提供区域多模型按权重预测（Node），提供半分布式水文模型搭建（Network）

## 这个包的使用风格应该是什么样子的？

`flux -> Element -> Node -> Network`

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

  model = Unit(
    name=name,
    surface=surf_ele,
    soil=soil_eles,
    route=route_ele,
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

## Unit的特性

- unit将分为三个基础层:surface,soil,route
- surface层可以对应于,深度学习模型的input层,其中可能就是简单的信息传输,也可能会存在融雪模块需要ode计算
- soil层可以包含多个element层,通常土壤的element都会认为时非线性水库,通常需要ode模块求解
- route层由于lagflux的限制,无法参与其他element进行联合计算,同样部分route层也会涉及ode计算

### Surface层特性

- surface层接受气象要素的输入,包括降雨,潜在蒸发,气温,日照时常等
- surface层的输出,全部统一名称为infiltration
- surface层还可以考虑一些不透水的情况,直接产出地表径流

### Soil层特性

- 受水箱模型的启发,soil层设计为可以包括多个element,设定为上下层等土壤分层
- soil层第一层的输入需要与surface层对接,所以接受变量必须为infiltration
- soil层会根据层数得出不同数量的出流量,每个层基本对应一个中间计算模块
- 
- 水文模型的核心计算模块，通常代表模型的基本输入
- 模型输出为flow，包括：
  1. 单一flow，如Exphydro
  2. Fastflow和Slowflow, 如GR4J
  3. Surfaceflow, Interflow, Baseflow三层
- 模型输入需要统一，比如Infiltration，Pet
