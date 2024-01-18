
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
