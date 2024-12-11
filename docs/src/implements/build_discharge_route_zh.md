# 构建河道汇流模型

## [discharge route model](https://gmd.copernicus.org/articles/14/7795/2021/)汇流模型介绍

河道汇流模型是描述水流在河道中运动过程的数学模型，主要用于计算河道中的流量变化。下面给出了一个简单的河道汇流模型的数学表达式:

```math
\begin{aligned}
\Delta S_{rf,n} &= Q_{rf,n-1} - Q_{rf,n} && (1) \\
Q_{rf,n} &= (S_{rf,n} + Q_{rf,n-1}) \cdot \frac{1}{\text{LAG}_{rf} + 1} && (2) \\
Q_{rf,0} &= \sum Q_{rf,up} && (3) \\
Q_{rf} &= Q_{rf,c_{max}} + (R_{sw} + R_{gw}) \cdot 0.001 \cdot A_{gc} && (4)
\end{aligned}
```

其中,$n$表示河段编号，取值范围为$\{1,...,c_{max}\}$;$\Delta S_{rf,n}$表示第n段河道的蓄水量变化;$Q_{rf,n}$表示第n段河道的出流量;$S_{rf,n}$表示第n段河道的蓄水量;$\text{LAG}_{rf}$ 表示河道汇流滞时参数;$Q_{rf,up}$表示上游河道的入流量;$R_{sw}$和$R_{gw}$分别表示地表径流和地下径流;$A_{gc}$表示栅格面积;$c_{max}$表示河道最大分段数.

公式(1)描述了河道蓄水量的变化，公式(2)表示河道出流量的计算，公式(3)给出了河道入流量的计算方法，公式(4)则是最终河道出流量的计算公式。
这种构建方式是根据每个hydrological response unit(HRU)的上游输入确定当前的HRU入流量,这个确定方式可以由两种路由函数来表示. 当前HRU的出流量则是以当前的河道蓄水量$S_{rf,n}$计算得到,并于HRU的R_{sw}和R_{gw}进行相加(经面积转换)得到最终的出流量$Q_{rf}$.

## HydroModels.jl框架下

### 计算公式的改造

HydroModels.jl框架中为这种基于河道蓄水状态计算时段出流量的方法以HydroRoute类型进行表示,通过构建ODE方程的方式对这个问题进行连续性求解.与原公式实现有所不同的是,,并认为在计算$Q_{rf,n}$时不存在瞬时变化,所以将公式(2)转换为如下式子:

```math
Q_{rf,n} &= S_{rf,n} \cdot \frac{1}{\text{LAG}_{rf} + 1} && (5) \\
```

式中,去掉了$Q_{rf,n-1}$项.因为认为,在持续性变化中,$S_{rf,n}$是一直处于变化状态的,所以在计算$Q_{rf,n}$时不需要考虑$Q_{rf,n-1}$的变化.此外于$S_{rf,n}$与$Q_{rf,n}$的量纲存在差异,同时如果将式中如果考虑$Q_{rf,n}$会额外增加一个状态变量,这个变量是一个瞬时值,也不适合作为状态变量来使用. 在HydroRoute构造中,其首要的就是针对出流计算公式进行表达,与HydroBucket一样,HydroRoute可以接受HydroFlux以此表示$Q_{rf,n}$的计算公式,如下所示:

```julia
rflux = HydroFlux([q, s_river] => [q_routed], [lag], exprs=[s_river / (1 + lag) + q])
```

即以Grid-Based和Vector-Based两种方式来表示. 在Route模块中是以