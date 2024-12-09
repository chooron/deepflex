function Base.show(io::IO, flux::AbstractHydroFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "HydroFlux(")
        print(io, "inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        print(io, ", params: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        print(io, ")")
    else
        println(io, "HydroFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        println(io, "  Expressions:")
        for (output, expr) in zip(flux.meta.outputs, flux.exprs)
            println(io, "    $output = $expr")
        end
    end

end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "StateFlux(")
        print(io, "inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        print(io, ", params: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        print(io, ", states: ", isempty(flux.meta.states) ? "nothing" : join(flux.meta.states, ", "))
        print(io, ")")
    else
        println(io, "StateFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        println(io, "  Parameters: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        println(io, "  States: ", isempty(flux.meta.states) ? "nothing" : join(flux.meta.states, ", "))
        println(io, "  Expressions:")
        println(io, "    $(flux.state) = $(flux.expr)")
    end
end

function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "NeuralFlux(")
        print(io, "inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        print(io, ", nn: ", join(flux.meta.nns, ", "))
        print(io, ")")
    else
        println(io, "NeuralFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        println(io, "  Neural Network: ", join(flux.meta.nns, ", "))
        println(io, "  Expression:")
        println(io, "    [$(join(flux.meta.outputs, ", "))] = $(flux.meta.nns[1])([$(join(flux.meta.inputs, ", "))])")
    end
end

function Base.show(io::IO, uh::AbstractHydrograph)
    compact = get(io, :compact, false)
    if compact
        print(io, "UnitHydroFlux(")
        print(io, "inputs: ", isempty(uh.meta.inputs) ? "nothing" : join(uh.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(uh.meta.outputs) ? "nothing" : join(uh.meta.outputs, ", "))
        print(io, ", params: ", isempty(uh.meta.params) ? "nothing" : join(uh.meta.params, ", "))
        print(io, ", uhfunc: ", nameof(typeof(uh.uhfunc).parameters[1]))
        print(io, ")")
    else
        println(io, "UnitHydroFlux:")
        println(io, "  Inputs: ", isempty(uh.meta.inputs) ? "nothing" : join(uh.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(uh.meta.outputs) ? "nothing" : join(uh.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(uh.meta.params) ? "nothing" : join(uh.meta.params, ", "))
        println(io, "  UnitFunction: ", nameof(typeof(uh.uhfunc).parameters[1]))
        println(io, "  SolveType: ", nameof(typeof(uh).parameters[end]))
    end
end

function Base.show(io::IO, ele::AbstractBucket)
    compact = get(io, :compact, false)
    if compact
        print(io, "Bucket(")
        print(io, "name: ", ele.meta.name)
        print(io, ", inputs: ", isempty(ele.meta.inputs) ? "nothing" : join(ele.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(ele.meta.outputs) ? "nothing" : join(ele.meta.outputs, ", "))
        print(io, ", params: ", isempty(ele.meta.params) ? "nothing" : join(ele.meta.params, ", "))
        print(io, ", states: ", isempty(ele.meta.states) ? "nothing" : join(ele.meta.states, ", "))
        print(io, ", funcs: ", length(ele.funcs), ", dfuncs: ", length(ele.dfuncs))
        print(io, ")")
    else
        println(io, "Bucket: ", ele.meta.name)
        println(io, "  Inputs: ", isempty(ele.meta.inputs) ? "nothing" : join(ele.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(ele.meta.outputs) ? "nothing" : join(ele.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(ele.meta.params) ? "nothing" : join(ele.meta.params, ", "))
        println(io, "  States: ", isempty(ele.meta.states) ? "nothing" : join(ele.meta.states, ", "))
        println(io, "  Flux functions: ", length(ele.funcs), " flux", length(ele.funcs) == 1 ? "" : "es")
        println(io, "  State derivative functions: ", length(ele.dfuncs), " flux", length(ele.dfuncs) == 1 ? "" : "es")
    end
end

function Base.show(io::IO, route::AbstractHydroRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroRoute(")
        print(io, "name: ", route.meta.name)
        print(io, ", inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        print(io, ", params: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        print(io, ", states: ", isempty(route.meta.states) ? "nothing" : join(route.meta.states, ", "))
        print(io, ", routefunc: ", typeof(route.rfunc))
        print(io, ")")
    else
        println(io, "HydroRoute: ", route.meta.name)
        println(io, "  Inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        println(io, "  States: ", isempty(route.meta.states) ? "nothing" : join(route.meta.states, ", "))
        println(io, "  RouteFunc: ", typeof(route.rfunc))
    end
end

function Base.show(io::IO, route::AbstractRapidRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "RapidRoute(")
        print(io, "name: ", route.meta.name)
        print(io, ", inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        print(io, ", params: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        print(io, ")")
    else
        println(io, "RapidRoute: ", route.meta.name)
        println(io, "  Inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
    end
end

function Base.show(io::IO, model::AbstractModel)
    fluxes_in_model = filter(x -> x isa AbstractFlux, model.components)
    buckets_in_model = filter(x -> x isa AbstractBucket, model.components)
    routes_in_model = filter(x -> x isa AbstractRoute, model.components)
    @assert(length(routes_in_model) <= 1, "Only one route is allowed in a model")
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroModel(")
        print(io, "name: ", model.meta.name)
        print(io, ", components: ", length(model.components))
        print(io, ")")
    else
        println(io, "HydroModel: ", model.meta.name)
        println(io, "  Components: ", join(map(c -> c.meta.name, model.components), ", "))
        println(io, "  Inputs: ", isempty(model.meta.inputs) ? "nothing" : join(model.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(model.meta.outputs) ? "nothing" : join(model.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(model.meta.params) ? "nothing" : join(model.meta.params, ", "))
        println(io, "  States: ", isempty(model.meta.states) ? "nothing" : join(model.meta.states, ", "))
        println(io, "  Components:")
        println(io, "    Fluxes: ", length(fluxes_in_model), " flux", length(fluxes_in_model) == 1 ? "" : "es")
        println(io, "    Buckets: ", length(buckets_in_model), " bucket", length(buckets_in_model) == 1 ? "" : "s")
        println(io, "    Route: ", isempty(routes_in_model) ? "nothing" : typeof(routes_in_model[1]))
    end
end