module DeepFlex

export ODEsElement

# 定义一个虚拟的Element类
abstract type BaseElement end
abstract type ParameterizedElement <: BaseElement end
abstract type StateElement <: BaseElement end
abstract type StateParameterizedElement <: BaseElement end
abstract type ODEsElement <: StateParameterizedElement end
abstract type DiscElement <: StateParameterizedElement end

end