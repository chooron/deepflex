# Element Methods
mutable struct ODEElement{T<:Number} <: AbstractElement
    id::String

    # parameters
    parameters::Dict{Symbol,T}
    param_names::Vector{Symbol}

    # states
    states::ComponentVector{T}
    init_states::ComponentVector{T}
    state_names::Vector{Symbol}
    multiplier::Dict{Symbol,ComponentVector{T}}

    # functions
    funcs::Vector{AbstractFunc}

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
end

function ODEElement(
    ; id::String,
    parameters::Dict{Symbol,T},
    init_states::ComponentVector{T},
    funcs::Vector{AbstractFunc},
    multiplier::Dict{Symbol,ComponentVector{T}}
)
    param_names = collect(keys(parameters))
    states = ComponentVector{T}(; Dict(k => [init_states[k]] for k in keys(init_states)))
    state_names = collect(keys(init_states))

    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        for nm in func.input_names
            if !(nm in input_names)
                push!(input_names, func.input_names)
            end
        end
        push!(output_names, func.output_name)
    end
    ODEElement(
        id,
        parameters,
        param_names,
        states,
        init_states,
        state_names,
        multiplier,
        funcs,
        input_names,
        output_names
    )
end

function get_id(ele::AbstractElement)::String
    return ele.id
end

function get_parameters(ele::AbstractElement; names::Vector{Symbol}=nothing)::Dict{Symbol,Any}
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function set_parameters!(ele::AbstractElement; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        ele.parameters[p.name] = p.value
        for func in ele.funcs
            set_parameters!(func, paraminfos=paraminfos)
        end
    end
end

function get_states(ele::AbstractElement; names::Vector{Symbol}=nothing)::Dict{Symbol,Vector{<:Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end

function set_states!(ele::Element; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        if p.name in ele.state_names
            setfield!(ele, p.name, p.value)
        end
    end
end

function get_init_states(ele::ODEElement)
    u_init_dict = Dict{Symbol,Number}()
    for sn in ele.state_names
        u_init_dict[sn] = getproperty(ele, sn)
    end
    return ComponentVector(u_init_dict)
end

function solve_prob(ele::ODEElement; input::ComponentVector{T}) where {T<:Number}
    dt = 1
    xs = 1:dt:length(input[first(keys(input))])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp = Dict(k => linear_interpolation(xs, input[k]) for k in keys(input))

    function ode_func!(du, u, p, t)
        # interpolate value by fitted functions
        tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
        tmp_du = get_du(ele, state=u, input=tmp_input)
        # return du
        du = ComponentVector(du; tmp_du...)
    end
    s_init = get_init_states(ele)
    prob = ODEProblem(ode_func!, s_init, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=dt)
    solved_u = sol.u
    solved_u_matrix = hcat(solved_u...)
    solved_u = ComponentVector(; Dict(nm => solved_u_matrix[idx, :] for (idx, nm) in enumerate(keys(solved_u[1])))...)
    return solved_u
end

function get_fluxes(ele::ODEElement; state::ComponentVector{T}, input::ComponentVector{T})
    fluxes = ComponentVector(input; state...)
    for func in ele.funcs
        temp_flux = get_output(func; input=fluxes)
        fluxes = ComponentVector(fluxes; temp_flux...)
    end
    return fluxes
end

function get_du(ele::ODEElement; state::ComponentVector{T}, input::ComponentVector{T})
    fluxes = get_fluxes(ele, state=state, input=input)
    du = ComponentVector(; Dict(k => [fluxes[pk] .* ele.weights[pk] for pk in keys(weight)] for (sk, weight) in ele.state_weights)...)
    return du
end

function get_output(ele::Element; input::ComponentVector{T}) where {T<:Number}
    state = solve_prob(ele, input=input)
    fluxes = get_fluxes(ele, state=state, input=input)
    return ComponentVector(; Dict(nm => fluxes[nm] for nm in ele.output_names))
end

# function get_output(ele::LagElement; input::Dict{Symbol,Vector{T}}) where {T<:Number}
#     weight = build_weight(ele)
#     states = solve_lag(weight, ele.lag_state, input)
#     # Get the new lag value to restart
#     final_states = states[end, :, :]
#     final_states[:, 1:end-1] = final_states[:, 2:end]
#     final_states[:, end] = 0
#     [states[:, i, 0] for i in 1:1:length(input[first(keys(input))])]
# end

# function solve_lag(weight::Vector{T}, lag_state::Vector{T}, input::Dict{Symbol,Vector{T}}) where {T<:Number}
#     max_length = max([len(w) for w in weight])

#     # output = np.zeros((len(input[0]), len(weight), max_length))  # num_ts, num_fluxes, len_lag

#     # for flux_num, (w, ls, i) in enumerate(zip(weight, lag_state, input)):
#     #     for ts in range(len(input[0])):
#     #         updated_state = ls + i[ts] * w
#     #         output[ts, flux_num, :len(w)] = updated_state[:]
#     #         ls = np.append(updated_state[1:], 0)
# end


function NeuralodeElement(in_dims::Int, out_dims::Int; node_layer)
    input_layer = Dense(in_dims, model_no_ode.layers[1].in_dims)
    output_layer = Dense(model_no_ode.layers[end].out_dims, out_dims)
    node_layer = node_layer
end