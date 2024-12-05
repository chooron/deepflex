#=using Main.HydroModels=#
:((var"##arg#2523971367920833662", var"##arg#7360549935612702713",
    var"##arg#10678911477305669588", var"##arg#5964805160111424296", var"##arg#14435632484259125414") -> begin
    #= D:\Julia\Julia-1.10.4\packages\packages\SymbolicUtils\EGhOJ\src\code.jl:373 =#
    #= D:\Julia\Julia-1.10.4\packages\packages\SymbolicUtils\EGhOJ\src\code.jl:374 =#
    #= D:\Julia\Julia-1.10.4\packages\packages\SymbolicUtils\EGhOJ\src\code.jl:375 =#
    begin
        prcp = var"##arg#2523971367920833662"[1]
        temp = var"##arg#2523971367920833662"[2]
        snowpack = var"##arg#7360549935612702713"[1]
        Tmin = var"##arg#10678911477305669588"[1]
        Tmax = var"##arg#10678911477305669588"[2]
        Df = var"##arg#10678911477305669588"[3]
        t = var"##arg#14435632484259125414"[1]
        begin
            snowfall = (*)((*)(0.5, prcp), (+)(1.0, (tanh)((*)(5.0, (+)(Tmin, (*)(-1, temp))))))
            rainfall = (*)((*)(0.5, prcp), (+)(1.0, (tanh)((*)(5.0, (+)((*)(-1, Tmin), temp)))))
            melt = (*)((*)(0.5, (+)(1.0, (tanh)((*)(5.0, (+)((*)(-1, Tmax), temp))))),
                (min)(snowpack, (*)(Df, (+)((*)(-1, Tmax), temp))))
            # used for ODE problem solving
            begin
                #= D:\Julia\Julia-1.10.4\packages\packages\SymbolicUtils\EGhOJ\src\code.jl:468 =#
                (SymbolicUtils.Code.create_array)(Vector, nothing, Val{1}(), Val{(1,)}(), (+)((*)(-1, melt), snowfall))
            end
            # used for flux calculation
            begin
                #= D:\Julia\Julia-1.10.4\packages\packages\SymbolicUtils\EGhOJ\src\code.jl:468 =#
                (SymbolicUtils.Code.create_array)(Vector, nothing, Val{1}(), Val{(3,)}(), snowfall, rainfall, melt)
            end
        end
    end
end)



RuntimeGeneratedFunction()
