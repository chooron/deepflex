function Base.show(io::IO, flux::AbstractHydroFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "HydroFlux(")
        print(io, "inputs: ", isempty(flux.meta.input) ? "nothing" : join(flux.meta.input, ", "))
        print(io, ", outputs: ", isempty(flux.meta.output) ? "nothing" : join(flux.meta.output, ", "))
        print(io, ", params: ", isempty(flux.meta.param) ? "nothing" : join(flux.meta.param, ", "))
        print(io, ")")
    else
        println(io, "HydroFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.input) ? "nothing" : join(flux.meta.input, ", "))
        println(io, "  Outputs: ", isempty(flux.meta.output) ? "nothing" : join(flux.meta.output, ", "))
        println(io, "  Parameters: ", isempty(flux.meta.param) ? "nothing" : join(flux.meta.param, ", "))
        println(io, "  Expressions:")
        for (output, expr) in zip(flux.meta.output, flux.exprs)
            println(io, "    $output = $expr")
        end
    end

end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "StateFlux(")
        print(io, "inputs: ", isempty(flux.meta.input) ? "nothing" : join(flux.meta.input, ", "))
        print(io, ", params: ", isempty(flux.meta.param) ? "nothing" : join(flux.meta.param, ", "))
        print(io, ", states: ", isempty(flux.meta.state) ? "nothing" : join(flux.meta.state, ", "))
        print(io, ")")
    else
        println(io, "StateFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.input) ? "nothing" : join(flux.meta.input, ", "))
        println(io, "  Parameters: ", isempty(flux.meta.param) ? "nothing" : join(flux.meta.param, ", "))
        println(io, "  States: ", isempty(flux.meta.state) ? "nothing" : join(flux.meta.state, ", "))
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

function Base.show(io::IO, flux::AbstractUnitHydroFlux)
    compact = get(io, :compact, false)
    if compact
        print(io, "UnitHydroFlux(")
        print(io, "inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        print(io, ", params: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        # print(io, ", uhfunc: ", nameof(flux.uhfunc))
        print(io, ")")
    else
        println(io, "UnitHydroFlux:")
        println(io, "  Inputs: ", isempty(flux.meta.inputs) ? "nothing" : join(flux.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(flux.meta.outputs) ? "nothing" : join(flux.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(flux.meta.params) ? "nothing" : join(flux.meta.params, ", "))
        # println(io, "  UnitHydrograph: ", nameof(flux.uhfunc))
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
        print(io, "GridRoute(")
        print(io, "name: ", route.meta.name)
        print(io, ", inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        print(io, ", outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        print(io, ", params: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        print(io, ", states: ", isempty(route.meta.states) ? "nothing" : join(route.meta.states, ", "))
        print(io, ", routefunc: ", typeof(route.rfunc))
        print(io, ", projtype: ", route.projtype)
        print(io, ", nodes: ", length(route.subareas))
        print(io, ")")
    else
        println(io, "GridRoute: ", route.meta.name)
        println(io, "  Inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        println(io, "  States: ", isempty(route.meta.states) ? "nothing" : join(route.meta.states, ", "))
        println(io, "  RouteFunc: StateRouteFlux(", typeof(route.rfunc).parameters[1], ")")
        println(io, "  ProjectionType: ", route.projtype)
        println(io, "  Nodes: ", length(route.subareas))
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
        print(io, ", routefunc: ", typeof(route.rfunc))
        print(io, ", nodes: ", size(route.adjacency)[1])
        print(io, ")")
    else
        println(io, "RapidRoute: ", route.meta.name)
        println(io, "  Inputs: ", isempty(route.meta.inputs) ? "nothing" : join(route.meta.inputs, ", "))
        println(io, "  Outputs: ", isempty(route.meta.outputs) ? "nothing" : join(route.meta.outputs, ", "))
        println(io, "  Parameters: ", isempty(route.meta.params) ? "nothing" : join(route.meta.params, ", "))
        println(io, "  RouteFunc: DynamicRouteFlux(", typeof(route.rfunc).parameters[1], ")")
        println(io, "  Nodes: ", size(route.adjacency)[1])
    end
end

function Base.show(io::IO, model::AbstractModel)
    fluxes_in_model = filter(x -> x isa AbstractFlux, model.components)
    buckets_in_model = filter(x -> x isa AbstractBucket, model.components)
    routes_in_model = filter(x -> x isa AbstractRoute, model.components)
    @assert(length(routes_in_model) <= 1, "Only one route is allowed in a model")
    model_type = :LumpedModel
    if length(routes_in_model) == 1
        if routes_in_model[1] isa AbstractGridRoute
            model_type = :GridModel
        elseif routes_in_model[1] isa AbstractVectorRoute
            model_type = :VectorModel
        elseif routes_in_model[1] isa AbstractSumRoute
            model_type = :MultiNodesModel
        else
            @error "Unknown route type: $(typeof(routes_in_model[1]))"
        end
    end

    compact = get(io, :compact, false)
    if compact
        print(io, "HydroModel(")
        print(io, "name: ", model.meta.name)
        print(io, ", type: ", model_type)
        print(io, ", components: ", length(model.components))
        print(io, ")")
    else
        println(io, "HydroModel: ", model.meta.name)
        println(io, "  Type: ", model_type)
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