function check_input(component::AbstractComponent, input::AbstractArray{<:Number,2}, timeidx::AbstractVector)
    @assert size(input, 1) == length(get_input_names(component)) "Input variables in component '$(component.meta.name)' do not match required dimensions. Expected $(length(get_input_names(component))) variables ($(get_input_names(component))), got $(size(input, 1)) variables"
    @assert size(input, 2) == length(timeidx) "Time steps in component '$(component.meta.name)' do not match required length. Expected $(length(timeidx)) steps, got $(size(input, 2)) steps"
end

function check_input(component::AbstractComponent, input::AbstractArray{<:Number,3}, timeidx::AbstractVector)
    @assert size(input, 1) == length(get_input_names(component)) "Input variables in component '$(component.meta.name)' do not match required dimensions. Expected $(length(get_input_names(component))) variables ($(get_input_names(component))), got $(size(input, 1)) variables"
    @assert size(input, 3) == length(timeidx) "Time steps in component '$(component.meta.name)' do not match required length. Expected $(length(timeidx)) steps, got $(size(input, 3)) steps"
end

function check_ptypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, ptypes::AbstractVector{Symbol})
    @assert length(ptypes) == size(input, 2) "Number of parameter types mismatch In $(get_name(component)). Expected $(size(input, 2)) parameter types, got $(length(ptypes))."
end

function check_stypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, stypes::AbstractVector{Symbol})
    @assert length(stypes) == size(input, 2) "Number of state types mismatch In $(get_name(component)). Expected $(size(input, 2)) state types, got $(length(stypes))."
end

function check_pas(component::AbstractComponent, pas::ComponentVector)
    check_parameters(component, pas)
    check_initstates(component, pas)
    check_nns(component, pas)
end

function check_pas(component::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{Symbol}, stypes::AbstractVector{Symbol})
    check_parameters(component, pas, ptypes)
    check_initstates(component, pas, stypes)
    check_nns(component, pas)
end

function check_parameters(component::AbstractComponent, pas::ComponentVector)
    param_names = get_param_names(component)
    cpt_name = get_name(component)
    for param_name in param_names
        @assert(param_name in keys(pas[:params]),
            "Parameter '$(param_name)' in component '$(cpt_name)' is required but not found in pas[:params]. Available parameters: $(keys(pas[:params]))"
        )
    end
end

function check_parameters(component::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{Symbol})
    param_names = get_param_names(component)
    cpt_name = get_name(component)
    for ptype in ptypes
        tmp_ptype_params_keys = keys(pas[:params][ptype])
        for param_name in param_names
            @assert(param_name in tmp_ptype_params_keys,
                "Parameter '$(param_name)' in component '$(cpt_name)' is required but not found in parameter type '$(ptype)'. Available parameters: $(tmp_ptype_params_keys)"
            )
        end
    end
end

function check_initstates(component::AbstractComponent, pas::ComponentVector)
    state_names = get_state_names(component)
    cpt_name = get_name(component)
    for state_name in state_names
        tmp_ptype_initstates_keys = keys(pas[:initstates])
        @assert(state_name in tmp_ptype_initstates_keys,
            "Initial state '$(state_name)' in component '$(cpt_name)' is required but not found in parameter type '$(ptype)'. Available states: $(tmp_ptype_initstates_keys)"
        )
    end
end

function check_initstates(component::AbstractComponent, pas::ComponentVector, stypes::AbstractVector{Symbol})
    state_names = get_state_names(component)
    cpt_name = get_name(component)
    for stype in stypes
        tmp_ptype_initstates_keys = keys(pas[:initstates][stype])
        for state_name in state_names
            @assert(state_name in tmp_ptype_initstates_keys,
                "Initial state '$(state_name)' in component '$(cpt_name)' is required but not found in state type '$(stype)'. Available states: $(tmp_ptype_initstates_keys)"
            )
        end
    end
end

function check_nns(component::AbstractComponent, pas::ComponentVector)
    nn_names = get_nn_names(component)
    cpt_name = get_name(component)
    if !isempty(nn_names)
        for nn_name in nn_names
            @assert nn_name in keys(pas[:nn]) "Neural network parameters '$(nn_name)' in component '$(cpt_name)' is required but not found in pas[:nn]. Available networks: $(keys(pas[:nn]))"
        end
    end
end