#=
* route 模块主要是涉及模块的汇流计算，就是如何将多个模块放到一起计算的
* routeflux是计算出每个节点汇流出去的结果，这个计算是依赖于当前节点的产流量（上游节点汇流加本节点汇流）
* 每个route之间的不同最本质的差别是获取上游汇流输入的不同
* 1. 在无上游汇入的话就是最基本的本节点的产流汇流
* 2. grid包括了上游输入，而上游输入方法是根据d8网络权重计算得到
* 3. vector包括了上游输入，而上游输入方法是根据有向图逐步计算得到
=#
"""
GridRoute是根据流向网格计算
"""
struct GridRoute <: AbstractRoute
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple
    """
    Hydrological flux functions
    """
    flux_func::Function
    """
    Hydrological ode functions
    """
    ode_func::Union{Nothing,Function}

    function GridRoute(
        name::Symbol;
        rfunc::Vector,
        funcs::Vector,
        dfuncs::Vector=StateFlux[],
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(funcs, dfuncs)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(vcat(funcs, dfuncs))
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(funcs)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, ode_func = build_ele_func(funcs, dfuncs, infos)

        return new(
            infos,
            flux_func,
            ode_func,
        )
    end
end