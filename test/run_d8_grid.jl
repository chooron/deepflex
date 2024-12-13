@testset "d8 grid routing" begin
    @testset "d8 grid routing (regular grid)" begin
        input = ones(9)
        positions = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
        flwdir = [1 4 8; 1 4 4; 1 1 2]
        result = HydroModels.grid_routing(input, positions, flwdir)
        target = [0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 2.0, 2.0]
        @test result == target
    end

    @testset "d8 grid routing (irregular grid)" begin
        input = ones(10)
        positions = [(1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 3), (4, 4)]
        flwdir = [0 4 4 0; 1 2 4 8; 0 1 4 0; 0 0 1 1]
        result = HydroModels.grid_routing(input, positions, flwdir)
        target = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 4.0, 1.0, 1.0]
        @test result == target
    end
end
