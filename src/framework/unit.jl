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

function pretrain!(unit::AbstractUnit; input::NamedTuple, train_config...)
    for ele in unit.elements
        pretrain!(ele; input=input, train_config...)
    end
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
    input::NamedTuple,
    parameters::NamedTuple,
    init_states::NamedTuple,
)
    # * This function is calculated element by element
    # initialize unit fluxes
    fluxes = input
    # traversal of the directed graph
    for ele in unit.elements
        if ele isa ODEElement
            solved_states = solve_prob(ele,
                input=fluxes, parameters=parameters,
                init_states=init_states[ele.state_names], solver=unit.solver[ele.name])
            tmp_fluxes = get_output(ele, input=fluxes, states=solved_states, parameters=parameters)
        else
            tmp_fluxes = get_output(ele, input=fluxes, parameters=parameters)
        end
        fluxes = merge(fluxes, tmp_fluxes)
    end
    fluxes
end