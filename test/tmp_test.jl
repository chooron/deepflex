using DataInterpolations
using ComponentArrays
using ModelingToolkit
using Symbolics

@variables t
function data_itp(t, time::AbstractVector, value::AbstractVector)
    itp = LinearInterpolation(value, time, extrapolate=true)
    itp(t)
end

input = ComponentVector(time=1:1:10, prcp=randn(10), temp=randn(10))
eqs = []
for nm in [:prcp, :temp]
    func_nm = Symbol(nm, "_itp")
    eval(Meta.parse("@variables $nm"))
    println(Meta.parse("$(func_nm)(t) = data_itp(t, time, input[$nm])").args)
    eval(Meta.parse("$(func_nm)(t) = data_itp(t, time, input[$nm])"))
    # eval(Expr(:=, Symbol("$func_nm(t)"), data_itp(t, time, input[nm])))
    eval(Meta.parse("@register_symbolic $func_nm(t)"))
    push!(eqs, eval(Meta.parse("$nm ~ $func_nm(t)")))
end

# @register_symbolic data_itp(t, time::AbstractVector, value::AbstractVector)
input_names = [:prcp, :temp]

# tmp_nm = :wind
# @variables $(input_names[1])
# function define_t()
#     vars = @variables t,v
#     vcat(vars, @variables prcp(t))
# end
# var = define_t()

# eqs = [var[3] ~ 1.0]
# test1 = ODESystem(eqs,name=:test1)

# eval(Meta.parse("snwr.$(:prcp)"))