# Coding Log

# 当前问题

- [X] 模型构建中需要调整fluxes的输入变量名称和输出变量名称
- [X] fluxes的function返回值过于固定，不够灵活
- [X] LAG Element无法参与整体ode的计算(可以参与计算但是兼容性很差,lagflux需要中间状态的缓冲计算)
- [X] SimpleFlux和HydroFlux无明显差异性
- [X] 构建多个function时，系统无法实现根据不同参数实现function的参数分配
- [ ] 计算成本相对于普通计算成本显著加大
  - [X] 原element计算产生的时间成本为2s, 直接构建一个简单的方程为38ms左右, 采用mtk.jl为40ms左右,另加mtk系统构建时间8ms
- [X] 各种ODE如何各自进行计算的话会加大插值产生的计算成本多余
  - [X] 采用mtk.jl对element公式进行整合
- [ ] lagflux如何改造成mtk.jl
  - [ ] lagflux改造成mtk的equations,在计算过程中其公式会不断发生改变
- [X] mtk.jl貌似只能支持一对一输入输出(已解决)
- [X] **由于component的参数存在多重嵌套，在参数优化的定义中存在问题**
- [X] 当前需要找出ODEProblem在用ForwardDifferetial求解时存在的问题，需要构建一个demo来重现这个问题，猜测这个问题应该是可调参数与不可调参数引起的问题
- [ ] 需要将水文的三种模型进行拆分，LumpedHydro.jl, SpatialHydro.jl, GridedHydro.jl

# 工作计划

- [X] routing function编写
- [X] 创建模型搭建基础类
- [X] 针对之前的模型进行ComponentArrays改造
- [X] 0实参构建模型
- [ ] 构建模型中，应该有构建合理性的校验和提示功能，主要的提示功能就是对模型各层计算数据是否存在缺失的功能，并对于缺失的模块提出增加建议
- [ ] 复现当前部分模型
- [ ] web端口构造
- [X] 在julia 1.10上完成部署
- [ ] StaticArrays或能够将性能进一步提升
- [ ] routing function的weight使用GuadGK.jl求解
- [X] 针对之前的模型进行ModelingToolkit改造
- [ ] **完善参数优化模块,包括模型参数优化,神经网络参数优化和混合参数优化**
- [ ] **提供自定义ODE求解,人为通过离散的方式求解,适应多数论文的计算,需要对比与DiscreteProblem之间的求解速度差距**
- [X] 将lag function嵌入至Node模块中
- [ ] Node中添加参数共享的设置
- [X] 将多个unit糅合到一块后应该如何表示参数，中间状态等参数，可以像lux.jl那样表示，比较清楚
- [ ] Node和Network的并行计算
- [ ] 搭建参数动态估计问题
- [ ] 如何将水文通量(Flux)转换成一系列类似于MTKstandarylibrary.jl那样的模块
- [X] input数据采用namedtuple类型,参数采用ComponentArray类型
- [X] superflexpy中unit是否具有存在意义, unit简单来说就是多个element的组合,可以考虑直接用elements list替代
- [X] 将input_names, output_names等信息通过函数调度，不作为模型存储的属性
- [X] Node将取消坡面汇流功能,将坡面汇流功能直接写在elements中
- [X] 参数信息的提取，可能需要进一步改进
- [X] 模型输入包括输入和参数,其中输入类型为NamedTuple,而参数为ComponentArray
- [X] 增加flux的图计算功能，保证即使flux顺序是混乱的也能够得到结果，并提供element增加出入通量的能力
- [ ] 构建node时提供构建信息

# 关键功能和实现技术

* [X] 基于julia多重分派的名称驱动代码风格，根据名称调用对应的flux
  * [X] 该项目是以命名驱动为核心，但是在一些情况中创建一些并没有实际含义的变量
* [X] Component基本计算实现, 涉及从某计算模块至整个流域汇流网络的计算
  * [X] 常微分方程求解(ModelingToolkit.jl, DifferentialEquations.jl)
  * [X] 网络拓扑计算(Graphs.jl)
* [X] Component基础的param_optimize参数优化功能实现
  * [X] 参数优化
* [ ] 参数动态模拟估计，Time-vary parameter estimation
* [X] 神经网络耦合物理公式计算的混合参数优化(包括普通参数和神经网络参数)

# 一些结论

* namedtuple类型合并相比componentArray要更加高效，而componentArray在兼容mtk时会更佳，在存储flux时采用namedtuple类型，而存储params采用componentVector类型
