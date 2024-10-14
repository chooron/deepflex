# Coding Log

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
- [X] routing function编写
- [X] 创建模型搭建基础类
- [X] 针对之前的模型进行ComponentArrays改造
- [X] 0实参构建模型
- [ ] 构建模型中，应该有构建合理性的校验和提示功能，主要的提示功能就是对模型各层计算数据是否存在缺失的功能，并对于缺失的模块提出增加建议
- [X] 复现当前部分模型
- [ ] web端口构造
- [X] 在julia 1.10上完成部署
- [ ] ~~StaticArrays或能够将性能进一步提升~~

- ~~ routing function的weight使用GuadGK.jl求解~~

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
- [ ] **~~根据macro提供一个自动生成模型计算的函数~~**
- [X] 根据计算网络结构迭代计算模型
- [X] LumpedHydro.jl中不考虑Node这个结构了，这个结构直接移至到SpatialHydro.jl
- [X] stateflux生成临时函数时存在问题
- [ ] sort_elements_by_topograph函数异常，或考虑不使用自动判断element计算顺序
- [X] 新增dPL-HBV, ENN, ~~PRNN~~
- [X] **~~NeuralFlux嵌入到dfunc无法生成耦合函数~~**
  - 当前输入变量只能是@varaibles (v(t))[1:4]这种类型，但这种类型或无法实现变量的替换
  - 考虑的方法是将nnflux前所有flux套入至nnflux中，但这种方式不行，因为nnflux前面可能还有nnflux
- [X] 我想让Flux的构建方式能够更有可读性，就是输入输出变量用键值对来连接

- ~~ 我记得当前在mtk框架下仍然难以通过AutoZygote的测试，这一块需要进一步完善~~

- [X] 非mtk框架下由于多次使用namedtuple，模型的计算性能还是不够好
- [ ] ~~记得本来采用StructArray，能够有效的避免反复计算带来的问题~~
- [ ] 自定义base.show
- [ ] ~~Zygote虽然不能用于mutable array, 但是可以通过chainrule执行自定义的rrule规则~~（不能对矩阵内部进行替换）
- [X] optimize需要提供多组数据训练的功能
- [ ] 当前optimization只针对于参数率定功能，后续可能会考虑
- [ ] 提供实时更新、添加、删除以及提示信息（包括当前element的输入输出）
- [X] 当前lagflux仅起到了信息记录作用，或可以删除=>改成routeflux可以用于各种汇流计算
- [X] neuralflux的参数或需要与其他参数独立出来，在分布式计算中不能对每个单元格都分配一个神经网络参数，故一般是一个统一的神经网络，所以参数类型为（ps=..., st=...., nn=...,）
- [ ] ~~模型输入的pas，三个主要键名：ps，st，nn~~
- [ ] 参数输入校验工作
- [ ] Route 类型的构建
  - [X] 马斯京跟
  - [X] 单位线
  - [X] hydrodischarge
- [ ] 使用macro构建simpleflux， @simpleflux var => expr, @lagflux var=> (flux, unithydro, lagtime), @stateflux var => expr, @neuralflux var => (input, nn) ，参考代码(temp/mtk_marco.jl)：
- [X] 当前DataInterpolations.jl在v5.0.0版本才兼容，其他版本存在问题
- [ ] DiscreteProblem似乎无法通过梯度优化
- [X] 构建了flux的计算匿名函数与state一致，**这时候就需要额外构建针对lag flux的element了**
- [X] 中间计算转为matrix合并的方式存储信息,效率显著提高
- [ ] 值得注意的是route flux可以分为两种，第一种类似于单位线这种，是在所有数据输入后才能计算出最终结果，第二种可以在一个时段后就可以得到汇流结果，比如马斯京根和cascade
- [X] 发现route这个过程可以理解成q_out转换为new_q_in的一个过程，不同就在于这个转换过程的不同
- [X] 完成了routeflux的核心构建，实现route模块与routeflux模块的拆分
- [X] 将routeflux细化成routeflux和unithydroflux两种类型
- [ ] vectorflux求解方式好像更倾向于discrete求解
- [ ] ~~将径流深(mm)转换为流量之后(m3/s),好像都不适用于continous solve, 这个需要进一步明确~~
- [ ] 运行需要额外的需求,比如:分段运行(保存上一次运行的中间状态); 持续运行(更新每次计算状态)
- [ ] 需要添加log信息和运行状态展示,用于前端展示
- [X] 马斯京根算法的连续性方程求解
- [ ] 插值方法这类的时变函数,可以设计成一个TimeVaryingFlux, 这样就可以直接参与连续性方程的求解了
- [ ] 设置基于TimeVaryingFlux的插值Flux的可行性,由于对于不同模块,需要重新构建插值函数了,所以插值函数这块我认为是没必要构建的

## 关键功能和实现技术

* [X] 基于julia多重分派的名称驱动代码风格，根据名称调用对应的flux
  * [X] 该项目是以命名驱动为核心，但是在一些情况中创建一些并没有实际含义的变量
* [X] Component基本计算实现, 涉及从某计算模块至整个流域汇流网络的计算
  * [X] 常微分方程求解(ModelingToolkit.jl, DifferentialEquations.jl)
  * [X] 网络拓扑计算(Graphs.jl)
* [X] Component基础的param_optimize参数优化功能实现
  * [X] 参数优化
* [ ] 参数动态模拟估计，Time-vary parameter estimation
* [X] 神经网络耦合物理公式计算的混合参数优化(包括普通参数和神经网络参数)


## 关于汇流计算的一些思考

vector-based河网汇流计算中, 设定了q_in, q_out, q_gen, 以及riv_st四种变量, 

在计算中:
- 设定上游流域为i,下游为j, q_out(i,t-1)应该是作为下游入流量q_in(j,t),
- 其中q_out的计算应该包括两个部分,首先是该流域的入流的汇流计算结果,即:
- q_out(i,t) = route_func(q_in(i,t), riv_st(i,t)) + q_gen(i,t)
- q_out可以根据流域入流和河床状态计算出流量,此外还应当包括流域本身的产流量.

在平衡计算中:
- q_out(i,t-1)就可以作为q_in(j,t)来替换q_in(j,t-1)
- riv_st(i,t)=riv_st(i,t-1)+q_in(i,t)-q_out(i,t), 注q_gen没有显式参与计算,但其包括在下游j的q_in中参与计算的
