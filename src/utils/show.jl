function Base.show(io::IO, flux::AbstractSimpleFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "SimpleFlux(")
        print(io, "inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        print(io, ", outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        print(io, ", params: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        print(io, ")")
    else
        println(io, "SimpleFlux:")
        println(io, "  Inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        println(io, "  Outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        println(io, "  Parameters: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        println(io, "  Expressions:")
        for (output, expr) in zip(flux.infos.output, flux.exprs)
            println(io, "    $output = $expr")
        end
    end

end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "StateFlux(")
        print(io, "inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        print(io, ", params: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        print(io, ", states: ", isempty(flux.infos.state) ? "nothing" : join(flux.infos.state, ", "))
        print(io, ")")
    else
        println(io, "StateFlux:")
        println(io, "  Inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        println(io, "  Parameters: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        println(io, "  States: ", isempty(flux.infos.state) ? "nothing" : join(flux.infos.state, ", "))
        println(io, "  Expressions:")
        println(io, "    $(flux.state) = $(flux.expr)")
    end
end

function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "NeuralFlux(")
        print(io, "inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        print(io, ", outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        print(io, ", nn: ", join(flux.infos.nn, ", "))
        print(io, ")")
    else
        println(io, "NeuralFlux:")
        println(io, "  Inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        println(io, "  Outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        println(io, "  Neural Network: ", join(flux.nnparam, ", "))
        println(io, "  Expression:")
        println(io, "    [$(join(flux.outputs, ", "))] = $(flux.infos.nn[1])([$(join(flux.inputs, ", "))])")
    end
end

function Base.show(io::IO, flux::AbstractRouteFlux)
    compact = get(io, :compact, false)
    rtype = typeof(flux).parameters[1]
    if compact
        print(io, "RouteFlux(")
        print(io, "inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        print(io, ", outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        print(io, ", params: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        print(io, ", states: ", isempty(flux.infos.state) ? "nothing" : join(flux.infos.state, ", "))
        print(io, ", rtype: ", rtype)
        print(io, ")")
    else
        println(io, "RouteFlux:")
        println(io, "  Inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        println(io, "  Outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        println(io, "  Parameters: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        println(io, "  States: ", isempty(flux.infos.state) ? "nothing" : join(flux.infos.state, ", "))
        println(io, "  RouteType: ", rtype)
    end
end

function Base.show(io::IO, flux::AbstractUnitHydroFlux)
    compact = get(io, :compact, false)
    if compact
        print(io, "UnitHydroFlux(")
        print(io, "inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        print(io, ", outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        print(io, ", params: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        print(io, ", uhfunc: ", nameof(flux.uhfunc))
        print(io, ")")
    else
        println(io, "UnitHydroFlux:")
        println(io, "  Inputs: ", isempty(flux.infos.input) ? "nothing" : join(flux.infos.input, ", "))
        println(io, "  Outputs: ", isempty(flux.infos.output) ? "nothing" : join(flux.infos.output, ", "))
        println(io, "  Parameters: ", isempty(flux.infos.param) ? "nothing" : join(flux.infos.param, ", "))
        println(io, "  UnitHydrograph: ", nameof(flux.uhfunc))
    end
end

function Base.show(io::IO, ele::AbstractBucket)
    println(io, ele.infos)
end