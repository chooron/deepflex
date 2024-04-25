#* 测试macro构建的interpolations系统
using DataInterpolations
using ModelingToolkit
using MacroTools
using ModelingToolkit: t_nounits as t

function build_data_itp_sys(input, time, varinfo)
    eqs = Equation[]
    for (nm, ip) in pairs(input)
        func_nm = Symbol(nm, :_itp)
        tmp_var = varinfo[nm]
        ex = quote
            $func_nm(t) = LinearInterpolation($ip, $time)(t)
            @register_symbolic $func_nm(t)
            push!($eqs, $tmp_var ~ $(func_nm)(t))
        end
        eval(ex)
    end
    eqs
end

p_data = rand(100)
e_data = rand(100)
l_data = rand(100)
input = (p=p_data, e=e_data, l=l_data)
@variables p(t) e(t) l(t)

equations = build_data_itp_sys(input, 1:100, (p=p, e=e, l=l))