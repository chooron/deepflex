# todo 像Sciml那样用Unit{false}和Unit{true}来区别基于mtk和非mtk的
struct Unit{E} <: AbstractUnit where {E<:AbstractElement}
    #* 名称信息
    name_info::NameInfo
    #* 表层计算模块
    surf_layer::E
    #* 土壤层计算模块
    soil_layer::Vector{E}
    #* 坡地汇流层计算模块
    slope_layer::E
end

function Unit(name::Symbol;
    surf_layer::E,
    soil_layer::Vector{E},
    slope_layer::E,
) where {E<:AbstractElement}
    #* 先从element中获取基础信息主要是输出,输入,参数名称等
    input_names, output_names, param_names, state_names = get_element_infos(vcat([surf_layer], soil_layer, [slope_layer]))
    #* 将信息整合到一个类中，便于管理
    name_info = NameInfo(name, input_names, output_names, param_names, state_names)
    Unit(
        name_info,
        surf_layer,
        soil_layer,
        slope_layer
    )
end

function (unit::Unit{false})(
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    # * This function is calculated element by element
    fluxes = input
    for ele in vcat([unit.surf_layer], unit.soil_layer, [unit.slope_layer])
        fluxes = merge(fluxes, ele(fluxes, params, init_states))
    end
    fluxes
end

struct MTKUnit{E} <: AbstractComponent where {E<:AbstractElement}
    #* 名称信息
    name_info::NameInfo
    #* 表层计算模块
    surf_layer::E
    #* 土壤层计算模块
    soil_layer::Vector{E}
    #* 坡地汇流层计算模块
    slope_layer::E
    #* variable and parameters
    var_info::NamedTuple
    param_info::NamedTuple
    #* system
    conn_eqs::Vector
    sys_dict::Dict{HydroElement,ODESystem}
end

function MTKUnit(name::Symbol;
    surf_layer::E,
    soil_layer::Vector{E},
    slope_layer::E,
) where {E<:AbstractElement}
    input_names, output_names, param_names, state_names = get_element_infos(vcat([surf_layer], soil_layer, [slope_layer]))
    name_info = NameInfo(name, input_names, output_names, state_names, param_names)
    var_info, param_info = init_var_param(name_info)
    conn_eqs, sys_dict = init_system(elements, var_info, param_info)

    MTKUnit{E}(
        name_info,
        surf_layer,
        soil_layer,
        slope_layer,
        var_info,
        param_info,
        conn_eqs,
        sys_dict,
    )
end

function init_system(elements::Vector{HydroElement}, var_info, param_info)
    #* prepare mtk building
    #* 构建参数和变量
    sys_dict = Dict(ele => build_ele_system(ele.funcs, ele.dfuncs, var_info, param_info, name=Symbol(ele.name, :_base_sys)) for ele in elements)
    ele_var_dict = Dict(ele => vcat(ele.input_names, ele.output_names, ele.state_names) for ele in elements)
    eqs = []
    #* 这里是遍历两次这个elements，就是通过对比两两数组的共同flux名称，然后将名称相同的名称构建方程
    for tmp_ele1 in elements
        for tmp_ele2 in elements
            if tmp_ele1.name != tmp_ele2.name
                share_var_names = intersect(ele_var_dict[tmp_ele1], ele_var_dict[tmp_ele2])
                for nm in share_var_names
                    push!(eqs, eval(Expr(:call, :~, getproperty(sys_dict[tmp_ele1], nm), getproperty(sys_dict[tmp_ele2], nm))))
                end
            end
        end
    end
    eqs, sys_dict
end

function solve_prob(
    unit::MTKUnit;
    input::NamedTuple,
    params::NamedTuple,
    init_states::NamedTuple,
)
    #* setup data
    itp_sys = build_itp_system(input, unit.var_info)
    eqs = copy(unit.conn_eqs)
    for (ele, sys) in pairs(unit.sys_dict)
        for nm in keys(input)
            if nm in ele.input_names
                push!(eqs, eval(Expr(:call, :~, getproperty(sys, nm), getproperty(itp_sys, nm))))
            end
        end
    end
    simple_sys = structural_simplify(compose(ODESystem(eqs, t; name=Symbol(unit.name, :_sys)), itp_sys, collect(values(unit.sys_dict))...))
    x0 = Pair{Num,eltype(init_states)}[]
    p = Pair{Num,eltype(params)}[]
    for (ele, sys) in pairs(unit.sys_dict)
        #* setup init states
        for nm in ele.state_names
            push!(x0, eval(Expr(:call, :(=>), getproperty(sys, nm), init_states[nm])))
        end
        #* setup parameters
        for nm in ele.param_names
            push!(p, eval(Expr(:call, :(=>), getproperty(sys, nm), params[nm])))
        end
    end
    prob = ODEProblem{true,SciMLBase.FullSpecialize}(simple_sys, x0, (1.0, length(input[:time])), p)
    sol = solve(prob, Tsit5(), saveat=input[:time])
    solved_u = hcat(sol.u...)
    state_names = collect(keys(init_states))
    namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)])
end

function (unit::MTKUnit)(input::NamedTuple, params::NamedTuple, init_states::NamedTuple)
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