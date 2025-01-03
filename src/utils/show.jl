function Base.show(io::IO, flux::AbstractHydroFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "HydroFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", outputs=[", join(get_output_names(flux), ","), "]")
        print(io, ", params=[", join(get_param_names(flux), ","), "]")
        print(io, ")")
    else
        println(io, "HydroFlux{$(typeof(flux).parameters[1])}:")
        println(io, "   Inputs: [", join(get_input_names(flux), ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(flux), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(flux), ", "), "]")
        
        if !isempty(flux.exprs)
            println(io, "  Expressions:")
            for (output, expr) in zip(flux.meta.outputs, flux.exprs)
                println(io, "    $output = ", expr)
            end
        end
    end
end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "StateFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", states=[", join(get_state_names(flux), ","), "]")
        print(io, ", params=[", join(get_param_names(flux), ","), "]")
        print(io, ")")
    else
        println(io, "StateFlux{$(typeof(flux).parameters[1])}:")
        println(io, "   Inputs: [", join(get_input_names(flux), ", "), "]")
        println(io, "   States: [", join(get_state_names(flux), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(flux), ", "), "]")
        
        if !isempty(flux.exprs)
            println(io, "  Expressions:")
            for (state, expr) in zip(flux.meta.states, flux.exprs)
                println(io, "    $state = ", expr)
            end
        end
    end
end

function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "NeuralFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", outputs=[", join(get_output_names(flux), ","), "]")
        print(io, ", nn=[", join(get_nn_names(flux), ","), "]")
        print(io, ")")
    else
        println(io, "NeuralFlux{$(typeof(flux).parameters[1])}:")
        println(io, "   Inputs: [", join(get_input_names(flux), ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(flux), ", "), "]")
        println(io, "   NNs: [", join(get_nn_names(flux), ", "), "]")
        
        if !isempty(flux.exprs)
            println(io, "  Expressions:")
            println(io, "    $(flux.meta.outputs) = $(flux.exprs)")
        end
    end
end

function Base.show(io::IO, uh::AbstractHydrograph)
    compact = get(io, :compact, false)
    if compact
        print(io, "UnitHydroFlux(")
        print(io, "inputs=[", join(get_input_names(uh), ","), "]")
        print(io, ", outputs=[", join(get_output_names(uh), ","), "]")
        print(io, ", params=[", join(get_param_names(uh), ","), "]")
        print(io, ", uhfunc: ", typeof(uh.uhfunc).parameters[1])
        print(io, ")")
    else
        println(io, "UnitHydroFlux:")
        println(io, "   Inputs: [", join(get_input_names(uh), ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(uh), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(uh), ", "), "]")
        println(io, "   UnitFunction: ", typeof(uh.uhfunc).parameters[1])
        println(io, "   SolveType: ", typeof(uh).parameters[end])
    end
end

function Base.show(io::IO, ele::AbstractBucket)
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroBucket(")
        print(io, "inputs=[", join(get_input_names(ele), ","), "]")
        print(io, ", states=[", join(get_state_names(ele), ","), "]")
        print(io, ", outputs=[", join(get_output_names(ele), ","), "]")
        print(io, ", params=[", join(get_param_names(ele), ","), "]")
        print(io, ")")
    else
        println(io, "HydroBucket{$(typeof(ele).parameters[1])}:")
        println(io, "   Inputs: [", join(get_input_names(ele), ", "), "]")
        println(io, "   States: [", join(get_state_names(ele) , ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(ele), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(ele), ", "), "]")
    end
end

function Base.show(io::IO, route::AbstractHydroRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroRoute(")
        print(io, "inputs=[", join(get_input_names(route), ","), "]")
        print(io, ", states=[", join(get_state_names(route), ","), "]")
        print(io, ", outputs=[", join(get_output_names(route), ","), "]")
        print(io, ", params=[", join(get_param_names(route), ","), "]")
        print(io, ")")
    else
        println(io, "HydroBucket{$(typeof(route).parameters[1])}:")
        println(io, "   Inputs: [", join(get_input_names(route), ", "), "]")
        println(io, "   States: [", join(get_state_names(route) , ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(route), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(route), ", "), "]")
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
        print(io, "name: ", get_name(model))
        print(io, ", components: ", length(model.components))
        print(io, ")")
    else
        println(io, "HydroModel: ", get_name(model))
        println(io, "   Components: ", join(map(c -> get_name(c), model.components), ", "))
        println(io, "   Inputs: [", join(get_input_names(model), ", "), "]")
        println(io, "   States: [", join(get_state_names(model) , ", "), "]")
        println(io, "   Outputs: [", join(get_output_names(model), ", "), "]")
        println(io, "   Parameters: [", join(get_param_names(model), ", "), "]")
        println(io, "   Components:")
        println(io, "       Fluxes: ", length(fluxes_in_model), " flux", length(fluxes_in_model) == 1 ? "" : "es")
        println(io, "       Buckets: ", length(buckets_in_model), " bucket", length(buckets_in_model) == 1 ? "" : "s")
        println(io, "       Routes: ", length(routes_in_model), " route", length(routes_in_model) == 1 ? "" : "s")
    end
end