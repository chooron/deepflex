# Coding Log

# 当前问题

- [X] 模型构建中需要调整fluxes的输入变量名称和输出变量名称
- [X] fluxes的function返回值过于固定，不够灵活
- [X] LAG Element无法参与整体ode的计算(可以参与计算但是兼容性很差,lagflux需要中间状态的缓冲计算)
- [X] SimpleFlux和HydroFlux无明显差异性
- [X] 构建多个function时，系统无法实现根据不同参数实现function的参数分配
- [X] 计算成本相对于普通计算成本显著加大
  - [X] 原element计算产生的时间成本为2s, 直接构建一个简单的方程为38ms左右, 采用mtk.jl为40ms左右,另加mtk系统构建时间8ms
- [X] 各种ODE如何各自进行计算的话会加大插值产生的计算成本多余
  - [X] 采用mtk.jl对element公式进行整合
- [ ] ~~lagflux如何改造成mtk.jl~~
  - [ ] ~~lagflux改造成mtk的equations,在计算过程中其公式会不断发生改变~~
- [X] mtk.jl貌似只能支持一对一输入输出(已解决)
- [X] **由于component的参数存在多重嵌套，在参数优化的定义中存在问题**
- [X] 当前需要找出ODEProblem在用ForwardDifferetial求解时存在的问题，需要构建一个demo来重现这个问题，猜测这个问题应该是可调参数与不可调参数引起的问题
- [X] 需要将水文的三种模型进行拆分，LumpedHydro.jl, SpatialHydro.jl, ~~GridedHydro.jl~~

# 工作计划

- [X] routing function编写
- [X] 创建模型搭建基础类
- [X] 针对之前的模型进行ComponentArrays改造
- [X] 0实参构建模型
- [ ] 构建模型中，应该有构建合理性的校验和提示功能，主要的提示功能就是对模型各层计算数据是否存在缺失的功能，并对于缺失的模块提出增加建议
- [X] 复现当前部分模型
- [ ] web端口构造
- [X] 在julia 1.10上完成部署
- [ ] ~~StaticArrays或能够将性能进一步提升~~
- [ ] routing function的weight使用GuadGK.jl求解
- [X] 针对之前的模型进行ModelingToolkit改造
- [X] **完善参数优化模块,包括模型参数优化,神经网络参数优化和混合参数优化**
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
- [ ] HBV计算结果有问题
- [X] M50无法实现在mtk下计算，以及step=false下计算
- [X] 将输入数据修改为StructArray类型

* [ ] **~~根据macro提供一个自动生成模型计算的函数~~**

- [X] 根据计算网络结构迭代计算模型
- [X] LumpedHydro.jl中不考虑Node这个结构了，这个结构直接移至到SpatialHydro.jl
- [X] stateflux生成临时函数时存在问题
- [ ] sort_elements_by_topograph函数异常，或考虑不使用自动判断element计算顺序
- [X] 新增dPL-HBV, ENN, ~~PRNN~~
- [X] **NeuralFlux嵌入到dfunc无法生成耦合函数**
  - 当前输入变量只能是@varaibles (v(t))[1:4]这种类型，但这种类型或无法实现变量的替换
  - 考虑的方法是将nnflux前所有flux套入至nnflux中，但这种方式不行，因为nnflux前面可能还有nnflux

# 暑假工作计划

前面有杂事耽搁了，最近可以重新优化这个框架啦，计划改造项目如下：

- [X] 我想让Flux的构建方式能够更有可读性，就是输入输出变量用键值对来连接
- [ ] 我记得当前在mtk框架下仍然难以通过AutoZygote的测试，这一块需要进一步完善
- [X] 非mtk框架下由于多次使用namedtuple，模型的计算性能还是不够好
- [ ] ~~记得本来采用StructArray，能够有效的避免反复计算带来的问题~~
- [ ] 自定义base.show
- [ ] ~~Zygote虽然不能用于mutable array, 但是可以通过chainrule执行自定义的rrule规则~~（不能对矩阵内部进行替换）
- [X] optimize需要提供多组数据训练的功能
- [ ] 当前optimization只针对于参数率定功能，后续可能会考虑
- [ ] 提供实时更新、添加、删除以及提示信息（包括当前element的输入输出）
- [ ] 使用macro构建simpleflux， @simpleflux var => expr, @lagflux var=> (flux, unithydro, lagtime), @stateflux var => expr, @neuralflux var => (input, nn) ，参考代码如下：

```julia
is_variable(x::Num) = is_variable(x.val)
function is_variable(x)
    if x isa ModelingToolkit.SymbolicUtils.Symbolic
        if haskey(x.metadata, ModelingToolkit.Symbolics.VariableSource)
            src = x.metadata[ModelingToolkit.Symbolics.VariableSource]
            return first(src) == :variables
        end
    end
    return false
end

- [ ] 当前DataInterpolations.jl在v5.0.0版本才兼容，其他版本存在问题
- [ ] DiscreteProblem似乎无法通过梯度优化
- [X] 构建了flux的计算匿名函数与state一致，**这时候就需要额外构建针对lag flux的element了**
- [ ] 中间计算转为matrix合并的方式存储信息,效率显著提高

@parameters a1
@variables a2
@variables b

eq1 = b ~ a1 + a2^2 + 1
eq1.lhs
Equation
get_variables(eq1)
nn = Lux.Chain(layer_1 = Dense(2 => 3, relu), name=:ann)

```

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

# 优化模型的问题

**AutoZygote**

- No MTK (Solved):
  **ERROR: MethodError: no method matching length(::ChainRulesCore.ZeroTangent)**
- MTK:
  **ERROR: Compiling Tuple{Type{Dict}, Dict{Any, Any}}: try/catch is not supported.**
