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
    compact = get(io, :compact, false)
    if compact
        print(io, "Bucket(")
        print(io, "name: ", ele.infos.name)
        print(io, ", inputs: ", isempty(ele.infos.input) ? "nothing" : join(ele.infos.input, ", "))
        print(io, ", outputs: ", isempty(ele.infos.output) ? "nothing" : join(ele.infos.output, ", "))
        print(io, ", params: ", isempty(ele.infos.param) ? "nothing" : join(ele.infos.param, ", "))
        print(io, ", states: ", isempty(ele.infos.state) ? "nothing" : join(ele.infos.state, ", "))
        print(io, ", funcs: ", length(ele.funcs), ", dfuncs: ", length(ele.dfuncs))
        print(io, ")")
    else
        println(io, "Bucket: ", ele.infos.name)
        println(io, "  Inputs: ", isempty(ele.infos.input) ? "nothing" : join(ele.infos.input, ", "))
        println(io, "  Outputs: ", isempty(ele.infos.output) ? "nothing" : join(ele.infos.output, ", "))
        println(io, "  Parameters: ", isempty(ele.infos.param) ? "nothing" : join(ele.infos.param, ", "))
        println(io, "  States: ", isempty(ele.infos.state) ? "nothing" : join(ele.infos.state, ", "))
        println(io, "  Flux functions: ", length(ele.funcs), " flux", length(ele.funcs) == 1 ? "" : "es")
        println(io, "  State derivative functions: ", length(ele.dfuncs), " flux", length(ele.dfuncs) == 1 ? "" : "es")
    end
end

function Base.show(io::IO, route::AbstractGridRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "GridRoute(")
        print(io, "name: ", route.infos.name)
        print(io, ", inputs: ", isempty(route.infos.input) ? "nothing" : join(route.infos.input, ", "))
        print(io, ", outputs: ", isempty(route.infos.output) ? "nothing" : join(route.infos.output, ", "))
        print(io, ", params: ", isempty(route.infos.param) ? "nothing" : join(route.infos.param, ", "))
        print(io, ", states: ", isempty(route.infos.state) ? "nothing" : join(route.infos.state, ", "))
        print(io, ", routefunc: ", typeof(route.rfunc))
        print(io, ", flwdir: ", size(route.flwdir))
        print(io, ", nodes: ", length(route.positions))
        print(io, ")")
    else
        println(io, "GridRoute: ", route.infos.name)
        println(io, "  Inputs: ", isempty(route.infos.input) ? "nothing" : join(route.infos.input, ", "))
        println(io, "  Outputs: ", isempty(route.infos.output) ? "nothing" : join(route.infos.output, ", "))
        println(io, "  Parameters: ", isempty(route.infos.param) ? "nothing" : join(route.infos.param, ", "))
        println(io, "  States: ", isempty(route.infos.state) ? "nothing" : join(route.infos.state, ", "))
        println(io, "  RouteFunc: RouteFlux(", typeof(route.rfunc).parameters[1], ")")
        println(io, "  FlowDirection: ", " Matrix{", eltype(route.flwdir), "}", size(route.flwdir))
        println(io, "  Nodes: ", length(route.positions))
    end
end

function Base.show(io::IO, route::AbstractVectorRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "VectorRoute(")
        print(io, "name: ", route.infos.name)
        print(io, ", inputs: ", isempty(route.infos.input) ? "nothing" : join(route.infos.input, ", "))
        print(io, ", outputs: ", isempty(route.infos.output) ? "nothing" : join(route.infos.output, ", "))
        print(io, ", params: ", isempty(route.infos.param) ? "nothing" : join(route.infos.param, ", "))
        print(io, ", states: ", isempty(route.infos.state) ? "nothing" : join(route.infos.state, ", "))
        print(io, ", routefunc: ", typeof(route.rfunc))
        println(io, ", rivernetwork: ", route.network)
        print(io, ", nodes: ", size(route.adjacency)[1])
        print(io, ")")
    else
        println(io, "VectorRoute: ", route.infos.name)
        println(io, "  Inputs: ", isempty(route.infos.input) ? "nothing" : join(route.infos.input, ", "))
        println(io, "  Outputs: ", isempty(route.infos.output) ? "nothing" : join(route.infos.output, ", "))
        println(io, "  Parameters: ", isempty(route.infos.param) ? "nothing" : join(route.infos.param, ", "))
        println(io, "  States: ", isempty(route.infos.state) ? "nothing" : join(route.infos.state, ", "))
        println(io, "  RouteFunc: RouteFlux(", typeof(route.rfunc).parameters[1], ")")
        println(io, "  RiverNetwork: $(summary(route.network)){$(nv(route.network)), $(ne(route.network))}")
        println(io, "  Nodes: ", size(route.adjacency)[1])
    end
end

function Base.show(io::IO, model::AbstractModel)
    # todo 
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
        print(io, "name: ", model.infos.name)
        print(io, ", type: ", model_type)
        print(io, ", components: ", length(model.components))
        print(io, ")")
    else
        println(io, "HydroModel: ", model.infos.name)
        println(io, "  Type: ", model_type)
        println(io, "  Inputs: ", isempty(model.infos.input) ? "nothing" : join(model.infos.input, ", "))
        println(io, "  Outputs: ", isempty(model.infos.output) ? "nothing" : join(model.infos.output, ", "))
        println(io, "  Parameters: ", isempty(model.infos.param) ? "nothing" : join(model.infos.param, ", "))
        println(io, "  States: ", isempty(model.infos.state) ? "nothing" : join(model.infos.state, ", "))
        println(io, "  Components:")
        println(io, "    Fluxes: ", length(fluxes_in_model), " flux", length(fluxes_in_model) == 1 ? "" : "es")
        println(io, "    Buckets: ", length(buckets_in_model), " bucket", length(buckets_in_model) == 1 ? "" : "s")
        println(io, "    Route: ", isempty(routes_in_model) ? "nothing" : typeof(routes_in_model[1]))
    end
end