struct Unit{E,S} <: AbstractUnit where {E<:AbstractElement,S<:AbstractSolver}
    name::Symbol

    # attribute
    input_names::Vector{Symbol}
    output_names::Vector{Symbol}
    state_names::Vector{Symbol}
    param_names::Vector{Symbol}

    # model structure
    elements::Vector{E}

    solver::Dict{Symbol,S}
end

function build_unit(; name::Symbol, elements::Vector{E}) where {E<:AbstractElement}
    input_names, output_names, state_names, param_names = get_element_info(elements)
    solver = Dict(ele.name => DEFAULT_SOLVER for ele in elements if ele isa ODEElement)

    Unit(
        name,
        input_names,
        output_names,
        state_names,
        param_names,
        elements,
        solver,
    )
end

function set_solver!(unit::Unit, solver::AbstractSolver)
    solver = Dict(ele.nm => solver for ele in unit.elements if ele isa ODEElement)
    setfield!(unit, :solver, solver)
end

function set_solver!(unit::Unit, solver::Dict{Symbol,S}) where {S<:AbstractSolver}
    setfield!(unit, :solver, solver)
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
    # * This function is calculated element by element
    # initialize unit fluxes
    fluxes = input
    # traversal of the directed graph
    for ele in unit.elements
        if ele isa ODEElement
            solved_states = solve_prob(ele,
                input=fluxes, parameters=parameters,
                init_states=init_states[ele.state_names], solver=unit.solver[ele.name])
            filter(isa(LagFlux), ele.funcs) do func
                init!(func, parameters)
            end
            tmp_fluxes = get_output(ele, input=fluxes, states=solved_states, parameters=parameters)
        else
            tmp_fluxes = get_output(ele, input=fluxes, parameters=parameters)
        end
        fluxes = ComponentVector(fluxes; tmp_fluxes...)
    end
    fluxes
end