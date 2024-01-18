# Element Methods
# get_id
function get_id(ele::BaseElement)::String
    return ele.id
end

# copy element
function copy_element(ele::ParameterizedElement)
    p = ele.parameters
    return typeof(ele)(ele.id, p)
end

function deepcopy_element(ele::ParameterizedElement)
    p = deepcopy(ele.parameters)
    return typeof(ele)(ele.id, p)
end

function copy_element(ele::StateElement)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, states)
end

function deepcopy_element(ele::StateElement)
    return copy_element(ele)
end

function copy_element(ele::StateElement)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, states)
end

function deepcopy_element(ele::StateElement)
    return copy_element(ele)
end

function copy_element(ele::Union{DiscElement,StateParameterizedElement})
    p = ele.parameters
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states)
end

function deepcopy_element(ele::Union{DiscElement,StateParameterizedElement})
    p = deepcopy(ele.parameters)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states)
end

function copy_element(ele::ODEsElement)
    p = ele.parameters
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states, ele.solver)
end

function deepcopy_element(ele::ODEsElement)
    p = deepcopy(ele.parameters)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states, ele.solver)
end


function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{Symbol}=nothing)::Dict{Symbol,Any}
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function get_states(ele::Union{StateElement,StateParameterizedElement}; names::Vector{Symbol}=nothing)::Dict{Symbol,Vector{<:Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end


function solve_prob(ele::ODEsElement; input::Dict{Symbol,Vector{T}})::Matrix{T} where {T<:Number}
    dt = 1
    xs = 1:dt:length(input[first(keys(input))])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp = Dict{Symbol, Any}()
    for (key, value) in pairs(input)
        itp[key] = linear_interpolation(xs, value)
    end

    function func(u, p, t)
        # interpolate value by fitted functions
        tmp_input = Dict{Symbol, T}()
        for (key, value) in pairs(itp)
            tmp_input[key] = value(t)
        end
        # return dt
        get_du(ele, u, tmp_input)
    end
    prob = ODEProblem(func, ele.init_states, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=dt)
    solved_u = sol.u
    return hcat(solved_u...)
end


function get_output(ele::ODEsElement; input::Dict{Symbol,Vector{T}}, solved::Bool=true)::Dict{Symbol,Vector{T}} where {T<:Number}
    S::Matrix{T} = solve_prob(ele, input=input)
    fluxes = get_fluxes(ele, S, input)
    return fluxes
end