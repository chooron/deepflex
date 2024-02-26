# 当前问题
- [ ] 模型构建中需要调整fluxes的输入变量名称和输出变量名称
- [ ] fluxes的function返回值过于固定，不够灵活

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