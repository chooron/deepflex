# Element Methods
mutable struct Element{T} <: AbstractElement where {T<:Number}
    name::String

    # parameters
    parameters::ComponentVector{T}

    # states
    states::ComponentVector{T}
    init_states::ComponentVector{T}
    state_names::Set{Symbol}

    # functions
    funcs::Vector{AbstractFlux}

    # solve config
    get_du::Function
    solve_type::String

    # attribute
    input_names::Set{Symbol}
    output_names::Set{Symbol}
end

function build_element(
    ; name::String,
    parameters::ComponentVector{T},
    init_states::ComponentVector{T},
    funcs::Vector{F},
    get_du::Function,
    solve_type::String=ode_solve,
) where {F<:AbstractFlux,T<:Number}

    states = ComponentVector(; Dict(k => [init_states[k]] for k in keys(init_states))...)
    state_names = Set(keys(init_states))

    input_names = Set{Symbol}()
    output_names = Set{Symbol}()

    for func in funcs
        union!(input_names, func.input_names)
        union!(output_names, func.output_names)
    end
    return Element{T}(
        name,
        parameters,
        states,
        init_states,
        state_names,
        funcs,
        get_du,
        solve_type,
        input_names,
        output_names
    )
end

function get_id(ele::AbstractElement)::String
    return ele.id
end

function get_parameters(ele::AbstractElement; names::Vector{Symbol}=nothing)
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function set_parameters!(ele::AbstractElement; paraminfos::Vector{P}) where {P<:AbstractParamInfo}
    for p in paraminfos
        ele.parameters[p.name] = p.value
        for func in ele.funcs
            set_parameters!(func, paraminfos=paraminfos)
        end
    end
end

function get_states(ele::AbstractElement; state_names::Union{Set{Symbol},Nothing}=nothing)
    if isnothing(state_names)
        return ele.states
    else
        available_state_names = [nm for nm in state_names if nm in ele.state_names]
        return ele.states[available_state_names]
    end
end


function pretrain!(ele::AbstractElement; input::ComponentVector{T}, train_config...) where {T<:Number}
    for func in ele.funcs
        if isa(func, LuxNNFunc)
            pretrain!(func, input=input, train_config...)
        end
    end
end

function solve_prob!(ele::AbstractElement; input::ComponentVector{T}) where {T<:Number}
    dt = 1
    xs = 1:dt:length(input[first(keys(input))])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp = Dict(k => linear_interpolation(xs, input[k]) for k in keys(input))

    function ode_func!(du, u, p, t)
        # interpolate value by fitted functions
        tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
        tmp_fluxes = get_fluxes(ele, state=u, input=tmp_input)
        tmp_du = ele.get_du(tmp_fluxes, ele.parameters)
        # return du
        for k in keys(ele.init_states)
            du[k] = tmp_du[k]
        end
    end

    prob = ODEProblem(ode_func!, ComponentVector(ele.init_states), tspan)
    sol = solve(prob, BS3(), dt=1.0, saveat=xs, reltol=1e-3, abstol=1e-3, sensealg=ForwardDiffSensitivity())
    solved_u = sol.u
    solved_u_matrix = hcat(solved_u...)
    solved_u = ComponentVector(; Dict(nm => solved_u_matrix[idx, :] for (idx, nm) in enumerate(keys(solved_u[1])))...)
    ele.states = solved_u
end

function get_fluxes(ele::AbstractElement; state::ComponentVector, input::ComponentVector{T}) where {T<:Number}
    fluxes = ComponentVector(input; state...)
    for func in ele.funcs
        temp_flux = func(fluxes)
        fluxes = ComponentVector(fluxes; temp_flux...)
    end
    return fluxes
end

function get_output(ele::AbstractElement; input::ComponentVector{T}) where {T<:Number}
    solve_prob!(ele, input=input)
    fluxes = get_fluxes(ele, state=ele.states, input=input)
    return ComponentVector(; Dict(nm => fluxes[nm] for nm in ele.output_names)...)
end

# function solve_prob!(ele::DCTElement; input::ComponentVector{T}) where {T<:Number}
#     data_len = length(input[first(keys(input))])
#     tmp_state = ele.init_states
#     for idx in 1:data_len
#         tmp_input = ComponentVector(; Dict(k => input[k][idx] for k in keys(input))...)
#         tmp_fluxes = get_fluxes(ele, state=tmp_state, input=tmp_input)
#         tmp_du = ele.get_du(tmp_fluxes, ele.parameters)
#         for k in keys(tmp_state)
#             update_state = tmp_state[k] + tmp_du[k]
#             tmp_state[k] = update_state
#             push!(ele.states[k], update_state)
#         end
#         tmp_state = ComponentVector(; Dict(k => tmp_state[k] + tmp_du[k] for k in keys(tmp_state))...)
#     end
# end
