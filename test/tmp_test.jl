using ModelingToolkit
using DifferentialEquations

struct TestModel
    @mtkmodel FOL begin
        @parameters begin
            Tmax
            Df
        end
        begin
            diff = Differential(t)
        end
        @variables begin
            t
            SnowWater(t)
            Snow(t)
            Temp(t)
        end
        @equations begin
            Melt ~ melt(SnowWater, Temp, Tmax, Df)
            diff(SnowWater) ~ SnowWater - Melt
        end
    end
end
using Lux
Lux.Chain