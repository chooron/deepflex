# Element Methods
# get_id
function get_id(ele::BaseElement)::String
    return ele.id
end

function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{Symbol}=nothing)::Dict{Symbol,Any}
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function set_parameters!(ele::Union{ParameterizedElement,StateParameterizedElement}; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        if p.name in ele.param_names
            setfield!(ele, p.name, p.value)
        end
    end
end

function get_states(ele::Union{StateElement,StateParameterizedElement}; names::Vector{Symbol}=nothing)::Dict{Symbol,Vector{<:Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end

function set_states!(ele::Union{StateElement,StateParameterizedElement}; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for p in paraminfos
        if p.name in ele.state_names
            setfield!(ele, p.name, p.value)
        end
    end
end

function get_output(ele::E; input::ComponentVector{T}) where {E<:BaseElement,T<:Number}
    fluxes = get_fluxes(ele, input=input)
    return fluxes
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

function get_init_states(ele::ODEsElement)
    u_init_dict = Dict{Symbol,Number}()
    for sn in ele.state_names
        u_init_dict[sn] = getproperty(ele, sn)
    end
    return ComponentVector(u_init_dict)
end

function solve_prob(ele::ODEsElement; input::ComponentVector{T}) where {T<:Number}
    dt = 1
    xs = 1:dt:length(input[first(keys(input))])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp = Dict(k=>linear_interpolation(xs, input[k]) for k in keys(input))

    function ode_func!(du, u, p, t)
        # interpolate value by fitted functions
        tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
        # return dt
        tmp_du = get_du(ele, S=u, input=tmp_input)
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


function get_output(ele::ODEsElement; input::ComponentVector{T}) where {T<:Number}
    S = solve_prob(ele, input=input)
    fluxes = get_fluxes(ele, S=S, input=input)
    return fluxes
end

@kwdef mutable struct LuxElement <: StateParameterizedElement
    # 模型参数
    model
    parameters
    states
    device
end

function LuxElement(model; device=cpu_device(), seed=42)
    rng = MersenneTwister()
    Random.seed!(rng, seed)
    ps, st = Lux.setup(rng, model) .|> device
    LuxElement(model=model, parameters=ps, states=st, device=device)
end

function LinearNNElement(in_dims::Int, out_dims::Int, hidd_size::Int, device=cpu_device(), seed=42)
    model = Chain(Dense(in_dims, hidd_size, tanh), Dense(hidd_size, out_dims))
    return LuxElement(model, device=device, seed=seed)
end

function update_lux_element!(ele::LuxElement, tstate)
    ele.model = tstate.model
    ele.parameters = tstate.parameters
    ele.states = tstate.states
end

function get_output(ele::LuxElement; input::ComponentVector{T}) where {T<:Number}
    x = hcat([input[k] for k in keys(input)]...)'
    y_pred = cpu_device()(Lux.apply(ele.model, ele.device(x), ele.parameters, ele.states)[1])
    return y_pred
end

function NeuralodeElement(in_dims::Int, out_dims::Int; node_layer)
    input_layer = Dense(in_dims, model_no_ode.layers[1].in_dims)
    output_layer = Dense(model_no_ode.layers[end].out_dims, out_dims)
    node_layer = node_layer
end