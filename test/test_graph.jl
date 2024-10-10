using Graphs

network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)

relation_mat = sparse(
    [src for (src, dst) in edges(network)],
    [dst for (src, dst) in edges(network)],
    1,
    nv(network), nv(network)
)

adj_matrix = adjacency_matrix(network)'

input_vec = Matrix(ones(1, 9))

reshape(adj_matrix * input_vec[1, :], 1, 9)