# `DeepFlex.jl` Code Log


# 当前问题
- [x] 模型构建中需要调整fluxes的输入变量名称和输出变量名称
- [x] fluxes的function返回值过于固定，不够灵活
- [x] LAG Element无法参与整体ode的计算(可以参与计算但是兼容性很差,lagflux需要中间状态的缓冲计算)
- [x] SimpleFlux和HydroFlux无明显差异性
- [x] 构建多个function时，系统无法实现根据不同参数实现function的参数分配


# 工作计划

- [ ] routing function编写
- [x] 创建模型搭建基础类
- [x] 针对之前的模型进行ComponentArrays改造
- [x] 直接将element分为三个大类：ODE，Lag和Simple
- [ ] 创建模型参数动态估计类，构建delay ode problem(一些参数是通过历史观测资料动态预测得到，而非通过优化算法根据实测预测结果无法计算得到)
- [x] 0实参构建模型
- [ ] 构建模型中，应该有构建合理性的校验和提示功能，主要的提示功能就是对模型各层计算数据是否存在缺失的功能，并对于缺失的模块提出增加建议
- [ ] 复现当前部分模型
- [ ] web端口构造
- [x] 在julia 1.10上完成部署
- [ ] StaticArrays或能够将性能进一步提升
- [ ] routing function的weight使用GuadGK.jl求解
- [x] 针对之前的模型进行ModelingToolkit改造, 准备做一个非mtl版本和mtl版本的包
- [ ] **完善参数优化模块,包括模型参数优化,神经网络参数优化和混合参数优化**
- [ ] **提供自定义ODE求解,人为通过离散的方式求解,适应多数论文的计算,需要对比与DiscreteProblem之间的求解速度差距**