
# macro combine_funcs(funcs...)
#     quote
#         function combined_func($(Symbol.(funcs[1].input_nms)...))
#             result = $(esc(funcs[1].func))($(Symbol.(funcs[1].input_nms)...))
#             $(for f in funcs[2:end]
#                 quote
#                     result = $(esc(f.func))($(Symbol.(f.input_nms)...), result)
#                 end
#             end)
#             return result
#         end
#         combined_func
#     end
# end

struct Fn
    input_nms
    output_nms
    func
end

# define two func
func1 = Fn([:a, :b], [:c], (input::NamedTuple) -> i[:a] + i[:b])
func2 = Fn([:a, :c], [:d], (input::NamedTuple) -> i[:a] * i[:c])

# 使用宏组合函数并返回
combined_func = @combine_funcs func1 func2

function combined_funcv2(input::NamedTuple)

end

# 测试组合函数
result = combined_func((2, 3))
println(result)  # 输出 10


expr_dict = Dict(
    :c => :(a + b),
    :d => :(c * b),
    :e => :(d + c)
)

# 已知变量
known_vars = Dict(:c => :(a + b), :d => :(c * b))

# 替换表达式中的变量
# 替换表达式中的变量
for (key, expr) in pairs(expr_dict)
    for (var, replacement) in known_vars
        expr = Meta.parse(string(expr))  # 将表达式转换为 Expr 类型
        expr = replace(expr, var => replacement)
    end
    expr_dict[key] = expr
end

combined_expr = :e
for (key, expr) in reverse(collect(expr_dict))
    combined_expr = :($(expr.args[1]) => $expr, $combined_expr)
end