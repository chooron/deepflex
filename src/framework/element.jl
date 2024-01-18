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

function get_states(ele::Union{StateElement,StateParameterizedElement}; names::Vector{Symbol}=nothing)::Dict{Symbol,Vector{<:Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end

function get_output(ele::E; input::Dict{Symbol,Vector{T}})::Dict{Symbol,Vector{T}} where {E<:BaseElement,T<:Number}
    fluxes = get_fluxes(ele, input=input)
    return fluxes
end

function solve_prob(ele::ODEsElement; input::Dict{Symbol,Vector{T}})::Matrix{T} where {T<:Number}
    dt = 1
    xs = 1:dt:length(input[first(keys(input))])
    tspan = (xs[1], xs[end])

    # fit interpolation functions
    itp = Dict{Symbol,Any}()
    for (key, value) in pairs(input)
        itp[key] = linear_interpolation(xs, value)
    end

    function func(u, p, t)
        # interpolate value by fitted functions
        tmp_input = Dict{Symbol,T}()
        for (key, value) in pairs(itp)
            tmp_input[key] = value(t)
        end
        # return dt
        get_du(ele, S=u, input=tmp_input)
    end
    prob = ODEProblem(func, ele.init_states, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=dt)
    solved_u = sol.u
    return hcat(solved_u...)
end


function get_output(ele::ODEsElement; input::Dict{Symbol,Vector{T}})::Dict{Symbol,Vector{T}} where {T<:Number}
    S::Matrix{T} = solve_prob(ele, input=input)
    fluxes = get_fluxes(ele, S=S, input=input)
    return fluxes
end