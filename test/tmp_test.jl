using Zygote

function f1(x)
    a = [x]
    sum(a)
end

gradient(f1, 1.0)