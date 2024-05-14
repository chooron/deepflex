function build_compute_topology(fluxes::AbstractVector{<:AbstractFlux})
    # 构建函数之间的计算图
    var_names = unique(vcat(get_input_output_names(fluxes)...))
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    topology = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(topology, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    topology
end

function build_compute_topology(elements::AbstractVector{<:AbstractElement})
    # 构建函数之间的计算图
    var_names = unique(vcat(get_input_output_names(elements)..., get_state_names(elements)...))
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    topology = SimpleDiGraph(length(var_names))
    for ele in elements
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(ele)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(topology, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    topology
end