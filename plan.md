# 工作计划

- [ ] routing function 编写
- [ ] 创建模型搭建基础类
- [ ] 针对之前的模型进行ModelingToolkit改造
- [x] 针对之前的模型进行ComponentArrays改造
- [ ] 直接将element分为三个大类：ODE，Dis和Lag，其他的公式计算就放入fluxes中


# 未来计划

- [ ] 复现当前部分模型
- [ ] web端口构造
- [ ] 在julia 1.10上完成部署
- [ ] StaticArrays或能够将性能进一步提升

开始复现模型时，发现重用性较差，一些element的公式仍旧写的过于笼统，使所有公式都进行拆分，element只体现flux的加减运算