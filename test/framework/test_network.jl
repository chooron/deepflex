using NamedTupleTools

include("../../src/network.jl")
include("../../src/utils/graph.jl")

struct MuskingumReach
    name::Symbol
    attr::NamedTuple
    params::NamedTuple
    upstream::Symbol
    downstream::Symbol
end
function (reach::MuskingumReach)(input::AbstractVector, dt::Number=1.0)
    output = zeros(eltype(input), size(input))
    output[1] = input[1]
    k, x = reach.params.k, reach.params.x

    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))

    for i in 1:(length(input)-1)
        q_out = (c0 * input[i+1]) + (c1 * input[i]) + (c2 * output[i])
        output[i+1] = max(minimum(input), q_out)
    end
    output
end

trange = 100
node_names = [:S1, :S2, :A1, :A2, :S3, :S4, :B1, :A3, :A4, :E1]
nodes = [(name=name,) for name in node_names]
result = namedtuple(
    node_names,
    [ones(trange) .* 7, ones(trange) .* 20, ones(trange) .* 27,
        ones(trange) .* 27, ones(trange) .* 10, ones(trange) .* 5,
        ones(trange) .* 15, ones(trange) .* 42, ones(trange) .* 42, ones(trange) .* 42])
reaches = [
    MuskingumReach(:reach_1, NamedTuple(), (k=1.0, x=0.2), :S1, :A1),
    MuskingumReach(:reach_2, NamedTuple(), (k=1.0, x=0.2), :S2, :A1),
    MuskingumReach(:reach_3, NamedTuple(), (k=1.0, x=0.2), :A1, :A2),
    MuskingumReach(:reach_4, NamedTuple(), (k=1.0, x=0.2), :A2, :A3),
    MuskingumReach(:reach_5, NamedTuple(), (k=1.0, x=0.2), :S3, :B1),
    MuskingumReach(:reach_6, NamedTuple(), (k=1.0, x=0.2), :S4, :B1),
    MuskingumReach(:reach_7, NamedTuple(), (k=1.0, x=0.2), :B1, :A3),
    MuskingumReach(:reach_8, NamedTuple(), (k=1.0, x=0.2), :A3, :A4),
    MuskingumReach(:reach_9, NamedTuple(), (k=1.0, x=0.2), :A4, :E1),
]

network = RiverNetwork(name=:net, nodes=nodes, reaches=reaches)
output = river_routing(network, result, 1.0)

# [node_names[idx] for idx in topological_sort(network.topology)]

