# todo 像Sciml那样用Unit{false}和Unit{true}来区别基于mtk和非mtk的
struct HydroUnit{E} <: AbstractUnit where {E<:AbstractElement}
    #* 名称信息
    nameinfo::NameInfo
    #* 表层计算模块
    surf_layer::Vector{E}
    #* 土壤层计算模块
    soil_layer::Vector{E}
    #* 坡地汇流层计算模块
    slope_layer::Vector{E}
end

function HydroUnit(name::Symbol;
    elements::Vector{E}
) where {E<:AbstractElement}
    #* 先从element中获取基础信息主要是输出,输入,参数名称等
    input_names, output_names, param_names, state_names = get_element_infos(elements)
    #* 将信息整合到一个类中，便于管理
    nameinfo = NameInfo(name, input_names, output_names, state_names, param_names)
    #* 根据elements的名称将其分配至unit的不同层中
    surf_layer, soil_layer, slope_layer = HydroElement[], HydroElement[], HydroElement[]
    for ele in elements
        if occursin("_surf", string(ele.nameinfo.name))
            push!(surf_layer, ele)
        elseif occursin("_soil", string(ele.nameinfo.name))
            push!(soil_layer, ele)
        elseif occursin("_slope", string(ele.nameinfo.name))
            push!(slope_layer, ele)
        else
            @error "$(ele.nameinfo.name) 不符合规范"
        end
    end
    HydroUnit(
        nameinfo,
        surf_layer,
        soil_layer,
        slope_layer
    )
end

function (unit::HydroUnit)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    # * This function is calculated element by element
    fluxes = input
    for ele in vcat(unit.surf_layer, unit.soil_layer, unit.slope_layer)
        fluxes = merge(fluxes, ele(fluxes, params, init_states))
    end
    fluxes
end

function get_all_luxnnflux(unit::HydroUnit)
    luxnn_tuple = namedtuple()
    for ele in vcat(unit.surf_layer, unit.soil_layer, unit.slope_layer)
        merge!(luxnn_tuple, get_all_luxnnflux(ele))
    end
    luxnn_tuple
end