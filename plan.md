# 当前问题
- [x] 模型构建中需要调整fluxes的输入变量名称和输出变量名称
- [x] fluxes的function返回值过于固定，不够灵活
- [ ] LAG Element无法参与整体ode的计算
- [ ] SimpleFlux和HydroFlux无明显差异性

# 创新点
- 耦合神经网络模型，提供一个PINN的水文模型构建框架
- 继承superflex和MARRMoT各自的优缺点，构建新的水文模型构建框架，具体来说：
    - 沿用superflex的构建框架，即element，node，network等
    - 沿用MARRMoT的flux function构建思路，使不同计算模块实现分类管理
    - **保证各模块的随机组合的可行性，当前计算模块可粗略分为SnowWater, SoilWater和RoutingStore三种**

## SoilWater的一些共性
- 水文模型的核心计算模块，通常代表模型的基本输入
- 模型输出为flow，包括：
    - (1) 单一flow，如Exphydro
    - (2) Fastflow和Slowflow, 如GR4J
    - (3) Surfaceflow, Interflow, Baseflow三层
- 模型输入需要统一，比如Infiltration，Pet


# 工作计划

- [ ] routing function 编写
- [x] 创建模型搭建基础类
- [x] 针对之前的模型进行ComponentArrays改造
- [x] 直接将element分为三个大类：ODE，Lag和Simple


# 未来计划

- [ ] 复现当前部分模型
- [ ] web端口构造
- [x] 在julia 1.10上完成部署
- [ ] StaticArrays或能够将性能进一步提升
- [ ] routing function 的weight使用GuadGK.jl求解
- [ ] 针对之前的模型进行ModelingToolkit改造, 准备做一个非mtl版本和mtl版本的包

开始复现模型时，发现重用性较差，一些element的公式仍旧写的过于笼统，使所有公式都进行拆分，element只体现flux的加减运算
构建多个function时，系统无法实现根据不同参数实现function的参数分配