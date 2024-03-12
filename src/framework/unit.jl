struct Unit{E,S} <: AbstractUnit where {E<:AbstractElement,S<:AbstractSolver}
    name::String

    # attribute
    input_names::Set{Symbol}
    output_names::Set{Symbol}
    state_names::Set{Symbol}
    parameter_names::Set{Symbol}

    # model structure
    elements::Vector{E}

    solver::Dict{Symbol,S}
end

function build_unit(; name::String, elements::Vector{E}, solver::Union{S,Dict{Symbol,S}}) where {E<:AbstractElement,S<:AbstractSolver}
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    parameter_names = Vector{Symbol}()

    if !(solver isa Dict)
        solver = Dict(ele.nm => solver for ele in elements if elements isa ODEElement)
    end

    for ele in elements
        union!(input_names, ele.input_names)
        union!(output_names, ele.output_names)
        union!(parameter_names, ele.parameter_names)
        if elements isa ODEElement
            union!(state_names, ele.state_names)
        end
    end

    Unit(
        name,
        input_names,
        output_names,
        state_names,
        parameter_names,
        elements,
        solver,
    )
end

function pretrain!(unit::AbstractUnit; input::ComponentVector{T}, train_config...) where {T<:Number}
    for ele in unit.elements
        pretrain!(ele; input=input, train_config...)
    end
end

function get_states(unit::AbstractUnit; state_names::Union{Set{Symbol},Nothing}=nothing)
    states = ComponentVector()
    for ele in unit.elements
        if ele isa ODEElement
            states = ComponentVector(states; get_states(ele, state_names=state_names)...)
        end
    end
    return states
end

function set_solver!(unit::Unit, solver::AbstractSolver)
    for ele in unit.elements
        if ele isa ODEElement
            set_solver!(ele, solver)
        end
    end
end

function get_output(
    unit::Unit;
    input::ComponentVector{T},
    parameters::ComponentVector{T},
    init_states::ComponentVector{T},
) where {T<:Number}
    # 开始计算
    # * This function is calculated element by element
    # initialize unit fluxes
    fluxes = input
    # traversal of the directed graph
    for tmp_ele in unit.elements
        if tmp_ele isa ODEElement
            solved_states = solve_prob(tmp_ele,
                input=fluxes, parameters=parameters, init_states=init_states,
                time_idx=input[:time], solver=unit.ele[tmp_ele.name])
            tmp_fluxes = get_output(tmp_ele, input=fluxes, states=solved_states, parameters=parameters)
        else
            tmp_fluxes = get_output(tmp_ele, input=fluxes, parameters=parameters)
        end
        fluxes = ComponentVector(fluxes; tmp_fluxes...)
    end
    fluxes
end