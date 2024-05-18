"""
$(TYPEDEF)


# Fields
$(FIELDS)
# Example
```

```
"""
struct MetaTopology
    "计算图"
    digraph::SimpleDiGraph
    "节点名称"
    node_names::Vector{Symbol}
    "输出名称与计算函数的映射"
    node_maps::NamedTuple
end

"""
$(SIGNATURES)

Construct calculation graphs based on all common hydrological fluxes in hydrological elements
"""
function build_compute_topology(fluxes::AbstractVector{<:AbstractFlux})
    var_names = unique(vcat(get_input_output_names(fluxes)...))
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    func_ntp = namedtuple(
        vcat([get_output_names(flux) for flux in fluxes]...),
        vcat([repeat([flux], length(get_output_names(flux))) for flux in fluxes]...)
    )

    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end

    MetaTopology(
        digraph,
        var_names,
        func_ntp
    )
end

"""
$(SIGNATURES)

Construct a calculation graph based on all hydrological elements in the hydrological unit
"""
function build_compute_topology(elements::AbstractVector{<:AbstractElement})
    var_names = unique(vcat(get_input_output_names(elements)..., get_state_names(elements)...))
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    ele_output_and_state_names(ele) = vcat(get_output_names(ele), get_state_names(ele))
    elements_ntp = namedtuple(
        vcat([ele_output_and_state_names(ele) for ele in elements]...),
        vcat([repeat([ele], length(ele_output_and_state_names(ele))) for ele in elements]...)
    )

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

    MetaTopology(
        digraph,
        var_names,
        elements_ntp
    )
end