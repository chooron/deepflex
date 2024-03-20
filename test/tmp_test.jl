using DataInterpolations
using ComponentArrays
using ModelingToolkit
using Symbolics

function data_itp(t, time::AbstractVector, value::AbstractVector)
    itp = LinearInterpolation(value, time, extrapolate=true)
    itp(t)
end

@register_symbolic data_itp(t, time::AbstractVector, value::AbstractVector)

input = ComponentVector(time=1:1:10, prcp=randn(10), temp=randn(10))


input_names = [:prcp, :temp]
# tmp_nm = :wind
# @variables $(input_names[1])
function define_t()
    vars = @variables t,v
    vcat(vars, @variables prcp(t))
end
var = define_t()

eqs = [var[3] ~ 1.0]
test1 = ODESystem(eqs,name=:test1)

eval(Meta.parse("snwr.$(:prcp)"))