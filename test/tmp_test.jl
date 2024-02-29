function bas!(a, c)
    a[:c] = c
    @info a
end

inp = (a = [1, 3], b = [2.5],)
bas!(inp, [4, 7])