
"""
Construct calculation graphs based on all common hydrological fluxes in hydrological components
"""
function sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})
    input_names = reduce(union, get_input_names.(fluxes))
    output_names = reduce(union, get_output_names.(fluxes))
    input_names = setdiff(input_names, output_names)
    output_names = setdiff(output_names, input_names)

    #* 构建flux输出名称与实例的namedtuple
    fluxes_ntp = reduce(merge, map(fluxes) do flux
        tmp_output_names = get_output_names(flux)
        NamedTuple{Tuple(tmp_output_names)}(repeat([flux], length(tmp_output_names)))
    end)

    #* 构建flux的有向计算图
    var_names = vcat(input_names, output_names)
    var_names_ntp = NamedTuple{Tuple(var_names)}(1:length(var_names))
    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                println((ipnm => opnm))
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    #* 根据有向图排序得到fluxes的计算顺序
    sorted_fluxes = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_flux = fluxes_ntp[tmp_var_nm]
            if !(tmp_flux in sorted_fluxes)
                push!(sorted_fluxes, tmp_flux)
            end
        end
    end
    sorted_fluxes
end

"""
Construct a calculation graph based on all hydrological components in the hydrological unit
"""
function sort_components(components::AbstractVector{<:AbstractComponent})
    input_names, output_names, state_names = get_var_names(components)
    components_ntp = reduce(merge, map(components) do component
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        println((tmp_input_names, tmp_output_names, tmp_state_names))
        tmp_output_state_names = vcat(tmp_output_names, tmp_state_names)
        NamedTuple{Tuple(tmp_output_state_names)}(repeat([component], length(tmp_output_state_names)))
    end)
    var_names = reduce(union, [input_names, output_names, state_names])
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    digraph = SimpleDiGraph(length(var_names))
    for component in components
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    sorted_components = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_component = components_ntp[tmp_var_nm]
            if !(tmp_component in sorted_components)
                push!(sorted_components, tmp_component)
            end
        end
    end
    sorted_components
end