using SparseArrays
using Zygote


# 定义矩阵的尺寸
n_rows = 3
n_cols = 4

# 定义对角线的值
main_diagonal = ones(n_rows)  # 主对角线值
upper_diagonal = ones(n_rows)  # 次对角线值

# 构建稀疏矩阵
sparse_matrix = spdiagm(0 => main_diagonal, 1 => upper_diagonal, 2 => upper_diagonal)
sum(sparse_matrix, dims=1)

function ff2(x)
    sparse_matrix = spdiagm([0 => x, 1 => x, 2 => x]...)
    sum(sparse_matrix)
end

gradient(ff2, [1, 2, 1])

sum(spdiagm([0 => [1, 2, 3], 1 => [2, 1, 2]]...), dims=1)