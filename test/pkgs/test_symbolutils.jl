using SymbolicUtils
using Symbolics: unwrap
function func1(t)
    t^2
end

func2(t) = t^2
func2(t::Num) = SymbolicUtils.term(func2,Symbolics.unwrap(t))
@variables a
func2(a)