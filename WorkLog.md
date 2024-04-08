# `DeepFlex.jl` Coding Log

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
- [ ] mtk.jl貌似只能支持一对一输入输出

# 工作计划

- [ ] routing function编写
- [X] 创建模型搭建基础类
- [X] 针对之前的模型进行ComponentArrays改造
- [X] 直接将element分为三个大类：ODE，Lag和Simple
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
- [ ] 将多个unit糅合到一块后应该如何表示参数，中间状态等参数，可以像lux.jl那样表示，比较清楚
- [ ] Node和Network的并行计算
- [ ] 搭建参数动态估计问题
