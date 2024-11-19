using Graphs
using Plots

flwdir = [0 4 8; 1 4 4; 1 1 2]
positions = [[1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]] # [1, 1], 

d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]

d8_nn_pads = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]

function build_grid_digraph(flwdir::AbstractMatrix, positions::AbstractVector)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]
    nrows, ncols = size(flwdir)
    # Count non-zero cells to determine number of nodes
    net = DiGraph(length(positions))

    # Helper function to convert row,col to linear index
    rc2idx(r, c) = (r - 1) * ncols + c

    # For each cell in grid
    for (i,(r, c)) in enumerate(positions)
        # Get flow direction code
        code = flwdir[r, c]

        # Skip if no flow direction
        code == 0 && continue

        # Find which D8 direction this matches
        dir_idx = findfirst(x -> x == code, d8_codes)
        isnothing(dir_idx) && continue

        # Calculate target cell based on padding pattern
        pad = d8_nn_pads[dir_idx]
        target_r = r + pad[1]
        target_c = c + pad[2]

        # Add edge if target is within bounds
        if 1 <= target_r <= nrows && 1 <= target_c <= ncols
            # Find target position index 
            target_idx = findfirst(p -> p == [target_r,target_c], positions)
            # Only add edge if both positions exist
            if !isnothing(target_idx)
                add_edge!(net, i, target_idx)
            end
        end
    end
    return net
end

net = build_grid_digraph(flwdir, positions)
adjacency_matrix(net)'
# Check network connectivity
println("Number of nodes: ", nv(net))
println("Number of edges: ", ne(net))
println("Network connectivity: ", edges(net))
# Print node connections
for e in edges(net)
    println("Node $(e.src) connects to node $(e.dst)")
end

# Print node connections with grid coordinates
for e in edges(net)
    src_r = div(e.src - 1, size(flwdir, 2)) + 1
    src_c = mod(e.src - 1, size(flwdir, 2)) + 1
    dst_r = div(e.dst - 1, size(flwdir, 2)) + 1
    dst_c = mod(e.dst - 1, size(flwdir, 2)) + 1
    println("Grid cell ($(src_r),$(src_c)) flows to ($(dst_r),$(dst_c))")
end
