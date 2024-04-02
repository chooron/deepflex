
struct MTKUnit{E} <: AbstractComponent where {E<:AbstractElement}
    #* 名称信息
    nameinfo::NameInfo
    #* 表层计算模块
    surf_layer::E
    #* 土壤层计算模块
    soil_layer::Vector{E}
    #* 坡地汇流层计算模块
    slope_layer::E
    #* variable and parameters
    varinfo::NamedTuple
    paraminfo::NamedTuple
end

function MTKUnit(
    name::Symbol;
    elements::Vector{E},
) where {E<:AbstractElement}
    #* 先从element中获取基础信息主要是输出,输入,参数名称等
    input_names, output_names, param_names, state_names = get_element_infos(elements)
    #* 将信息整合到一个类中，便于管理
    nameinfo = NameInfo(name, input_names, output_names, param_names, state_names)
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
    varinfo, paraminfo = init_var_param(nameinfo)
    base_systems = Dict(nm => [build_ele_system(ele, varinfo=varinfo, paraminfo=paraminfo) for ele in layer]
                        for (nm, layer) in zip([:surf, :soil, :slope], [surf_layer, soil_layer, slope_layer]))
    conn_eqs = init_system(elements)

    MTKUnit{E}(
        nameinfo,
        surf_layer,
        soil_layer,
        slope_layer,
        varinfo,
        paraminfo,
        conn_eqs,
        sys_dict,
    )
end

function get_all_elements(unit::MTKUnit)
    vcat(unit.surf_layer, unit.soil_layer, unit.slope_layer)
end

function get_ele_connect_eqs(unit::MTKUnit)
    #* prepare mtk building
    #* 构建参数和变量
    eqs = []
    element_list = get_all_elements(unit)
    ele_var_dict = Dict(ele => vcat(ele.input_names, ele.output_names, ele.state_names) for ele in element_list)
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in element_list
        for tmp_ele2 in filter(ele -> ele.name != tmp_ele1.name, element_list)
            share_var_names = intersect(ele_var_dict[tmp_ele1], ele_var_dict[tmp_ele2])
            for nm in share_var_names
                push!(eqs, eval(Expr(:call, :~, getproperty(tmp_ele1.base_sys, nm), getproperty(tmp_ele1.base_sys, nm))))
            end
        end
    end
    eqs
end

function setup_data(
    unit::MTKUnit;
    input::NamedTuple,
    time::AbstractVector
)
    itp_sys = build_itp_system(input, time, unit.varinfo, name=unit.nameinfo.name)
    eqs = get_ele_connect_eqs(unit)
    element_list = get_all_elements(unit)
    for ele in element_list
        for nm in filter(nm -> nm in ele.input_names, keys(input))
            push!(eqs, eval(Expr(:call, :~, getproperty(ele.base_sys, nm), getproperty(itp_sys, nm))))
        end
    end
    structural_simplify(
        compose(ODESystem(eqs, t; name=Symbol(unit.name, :_sys)),
            itp_sys, [ele.base_sys for ele in element_list]...)
    )
end

function solve_prob(
    unit::MTKUnit;
    sys::ODESystem,
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    x0 = Pair{Num,eltype(init_states)}[]
    p = Pair{Num,eltype(params)}[]
    for ele in get_all_elements(unit)
        #* setup init states
        for nm in ele.nameinfo.state_names
            push!(x0, eval(Expr(:call, :(=>), getproperty(ele.base_sys, nm), init_states[nm])))
        end
        #* setup parameters
        for nm in ele.nameinfo.param_names
            push!(p, eval(Expr(:call, :(=>), getproperty(ele.base_sys, nm), params[nm])))
        end
    end
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(sys, x0, (1.0, length(input[:time])), p)
    sol = solve(prob, Tsit5(), saveat=input[:time])
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)])
end

function (unit::MTKUnit)(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple
)
    # //先判断在slope层中是否存在lagfunc,如果存在则需要将slope层单独进行计算
    #* 直接考虑slope层与其他两个层单独进行计算
    
    #* 先从通过unit构建前两层的综合ODE问题
    states = solve_prob(unit, input=input, params=params, init_states=init_states)
    #* 将求解的结果带入到方程里面进行求解
    fluxes = merge(input, states)
    for ele in vcat([unit.surf_layer], unit.soil_layer)
        for func in ele.funcs
            fluxes = merge(fluxes, func(fluxes, params))
        end
    end
    merge(fluxes, unit.slope_layer(fluxes, params))
end