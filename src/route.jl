#=
* route 模块主要是涉及模块的汇流计算，就是如何将多个模块放到一起计算的
* 分布式模型中是根据
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