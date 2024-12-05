using Zygote
using BenchmarkTools

function single_iteration(state, p)
    # Perform your calculation here
    # This is where you would use the input mat and create a new mat
    mat, history = state
    new_mat = mat .* p
    return new_mat, vcat(history, [new_mat])
end

mat = ones(3, 3)
final_mat, history = reduce((acc, _) -> single_iteration(acc, 2), 1:10, init=(mat, (mat,)))



