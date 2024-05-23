"""
$(SIGNATURES)

Construct calculation graphs based on all common hydrological fluxes in hydrological elements
"""
function sort_fluxes_by_topograph(fluxes::AbstractVector{<:AbstractFlux})
    input_names, output_names = get_input_output_names(fluxes)
    var_names = vcat(input_names, output_names)
    var_names_ntp = NamedTuple{Tuple(var_names)}(1:length(var_names))

    fluxes_output_names, fluxes_funcs = [], []
    for flux in fluxes
        temp_flux_output_names = get_output_names(flux)
        push!(fluxes_output_names, temp_flux_output_names)
        push!(fluxes_funcs, repeat([flux], length(temp_flux_output_names)))
    end
    func_ntp = NamedTuple{Tuple(vcat(fluxes_output_names...))}(vcat(fluxes_funcs...))

    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end

    sorted_fluxes = AbstractFlux[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_flux = func_ntp[tmp_var_nm]
            if !(tmp_flux in sorted_fluxes)
                push!(sorted_fluxes, tmp_flux)
            end
        end
    end
    sorted_fluxes
end

"""
$(SIGNATURES)

Construct a calculation graph based on all hydrological elements in the hydrological unit
"""
function sort_elements_by_topograph(elements::AbstractVector{<:AbstractElement})
    input_names, output_names, state_names = get_var_names(elements)
    var_names = vcat(input_names, output_names, state_names)
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))

    elements_output_names, elements_structs = [], []
    for ele in elements
        temp_ele_output_names = vcat(get_var_names(ele)[[2, 3]]...)
        push!(elements_output_names, temp_ele_output_names)
        push!(elements_structs, repeat([ele], length(temp_ele_output_names)))
    end
    elements_ntp = NamedTuple{Tuple(vcat(elements_output_names...))}(vcat(elements_structs...))

    digraph = SimpleDiGraph(length(var_names))
    for ele in elements
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(ele)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end

    sorted_elements = AbstractElement[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_ele = elements_ntp[tmp_var_nm]
            if !(tmp_ele in sorted_elements)
                push!(sorted_elements, tmp_ele)
            end
        end
    end
    sorted_elements
end