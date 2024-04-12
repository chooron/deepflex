include("../../src/reach.jl")
using BenchmarkTools

bottom_width, side_slope, mannings_n, slope, reach_length = 100, 100, 0.1, 0.001, 2000
reach = TrapezoidalReach(bottom_width, reach_length, slope, mannings_n, side_slope)
input = [1.0, 64.44444444, 124.44444444, 180.0,
    231.11111111, 277.77777778, 320.0, 357.77777778,
    391.11111111, 420.0, 444.44444444, 464.44444444,
    480.0, 491.11111111, 497.77777778, 500.0,
    497.77777778, 491.11111111, 480.0, 464.44444444,
    444.44444444, 420.0, 391.11111111, 357.77777778,
    320.0, 277.77777778, 231.11111111, 180.0,
    124.44444444, 64.44444444, 1.0, 1.0,]

c0, c1, c2 = calcu_muskingum_params(reach, 100, 0.1)
#! mtk based problem 相比 common problem计算效率更低，因此选用common problem用于计算
@btime output = reach(input, 0.1)