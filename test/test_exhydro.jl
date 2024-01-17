function my_function(; kwargs...)
    # kwargs 是一个 Dict 类型的变量，包含了所有传入的关键字参数
    for (key, value) in kwargs
        println("Key: $key, Value: $value")
    end
end

# 使用关键字参数调用函数