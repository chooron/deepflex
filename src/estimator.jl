"""
# 参数估计的输入类型是比较多样的,包括单节点输入,时段输入
# 这个类的面向的问题在于:
# 1. 初始状态的确定,一般来说在参数优化中可能会考虑初始状态的影响,
#     但是比如土壤初始含水量的估计是与其最大含水量挂钩的,所以一般来说是用初始百分比来确定的,这时候就需要estimator设置一个百分比的参数
# 2. 然后就是一些参数的确定,比如马斯京根算法k,这些是可以直接根据河道长度推算的,所以可以采用estimator设置输入参数为河长,依次估计k
# 3. 还有一些参数的估计是根据数据确定的,比如说前期降雨序列,通过这个输入到神经网络中可以估计初始土壤含水

# 总之这个estimator是为了避免在原模型上大幅度魔改,使其有规范性,且降低修改门槛,灵活性也更强
# estimator的实质是将计算结果替换参数或状态,而其余类型都是计算结果给出新的变量
"""
struct HydroEstimator <: AbstractEstimator
    "用于预测的方法，一般来说一个Flux就够了"
    eflux::AbstractFlux
    "待用于预测的params"
    params::Vector{Num}
    "MetaData"
    infos::NamedTuple

    function HydroEstimator(eflux::AbstractFlux, params::Vector{Num})
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(eflux)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(eflux)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(eflux)
        #* Setup the name information of the hydrobucket
        infos = (name=name, input=input_names, output=output_names, state=state_names, param=param_names, nn=nn_names)

        return new(eflux, params, infos)
    end
end

function (est::HydroEstimator)!(input::AbstractArray, params::ComponentVector)
    #* 计算后需要执行替换工作
    updated_params = est.eflux(input, params)
    #* 生成用于合并的ComponentVector
    return est.eflux(input, params)
end


