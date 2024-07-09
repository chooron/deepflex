using Zygote

function test_func(i::NamedTuple, p::NamedTuple)
    # tmp_i = (c=,)
    # i_list = [i[k] for k in (:a,:b)]
    i_list = collect(i)
    tmp_i = NamedTuple{(:c,)}([i_list[1] .+ i_list[2] .^ 2 .* p.a])
    i = merge(i, tmp_i)
    sum(i.c .+ i.b .+ i.a)
end

Zygote.gradient(test_func, (a=[2, 3], b=[3, 3]), (a=[2, 1],))

v = (a=[2, 3], b=[3, 3])
test_func(v, (a=[2, 1],))
using Zygote

function loss_adjoint(p)
    prediction = p .* rand(2, 100)
    prediction2 = [prediction[:, i] for i in axes(prediction, 2)]
    sum(sum.(prediction2))
end

Zygote.gradient(loss_adjoint, ones(2))