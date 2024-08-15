"""
$(TYPEDEF)
The basic hydrological calculation unit usually contains multiple hydrological calculation modules,
Each hydrological calculation module can be solved through step-by-step calculation and overall calculation.
It is usually used to simulate the vertical calculation process of a watershed or a calculation unit.

a basic hydrology unit must include:
    Surface water layer, typical Elements include snowfall module, interception module, evaporation module, infiltration module, etc.
    Soil, the water layer in the soil. Typical Elements include soil moisture module, runoff calculation module, etc.
    FreeWater layer, typical elements include groundwater, surface water, soil flow, etc.

# Fields
$(FIELDS)
# Example
```
using LumpedHydro
using LumpedHydro.ExpHydro: Surface, Soil, FreeWater

HydroUnit(
    name,
    elements=[Surface(name=name, mtk=mtk), Soil(name=name, mtk=mtk), FreeWater(name=name, mtk=mtk)],
    step=step,
)
```
"""
struct HydroUnit <: AbstractUnit
    "hydrological computation unit information"
    infos::NamedTuple
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "input idx for each components"
    input_idx::Vector

    function HydroUnit(name; components::Vector{<:AbstractComponent})
        #* 获取每个element的输出结果,然后与输入结果逐次拼接,获取每次输入的matrix的idx
        unit_input_names = unit_var_names = get_var_names(components)[1]
        input_idx = Vector[]
        for component in components
            tmp_input_idx = map(get_input_names(component)) do nm
                findfirst(varnm -> varnm == nm, unit_var_names)
            end
            #* 更新unit_var_names
            unit_var_names = reduce(vcat, [unit_var_names, get_state_names(component), get_output_names(component)])
            push!(input_idx, tmp_input_idx)
        end
        unit_infos = (name=name, input=unit_input_names, var=unit_var_names)
        new(
            unit_infos,
            components,
            input_idx,
        )
    end
end

# 求解并计算
function (unit::HydroUnit)(
    input::NamedTuple,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver(),
    output_type::Symbol=:namedtuple
)
    fluxes = reduce(hcat, [input[nm] for nm in get_input_names(unit)])'
    for (tmp_ele, idx) in zip(unit.components, unit.input_idx)
        tmp_fluxes = tmp_ele(fluxes[idx, :], pas, timeidx=timeidx, solver=solver)
        fluxes = cat(fluxes, tmp_fluxes, dims=1)
    end
    if output_type == :namedtuple
        return NamedTuple{Tuple(get_var_names(unit))}(eachrow(fluxes))
    elseif output_type == :array
        return fluxes
    end
end

#* 多输入构建大型方程求解并计算
function (unit::HydroUnit)(
    inputs::Vector{<:NamedTuple},
    pas::ComponentVector;
    timeidx::Vector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
    solver::AbstractSolver=ODESolver(),
    output_type::Symbol=:array
)
    fluxes = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, [input[nm] for nm in get_input_names(unit)]) for input in inputs])
    fluxes = permutedims(fluxes, (2, 3, 1))
    for (ele, idx) in zip(unit.components, unit.input_idx)
        fluxes_input = fluxes[idx, :, :]
        if !isnothing(get_ode_func(ele))
            #* Call the solve_prob method to solve the state of element at the specified timeidx
            solved_states = solve_multi_prob(
                ele, input=fluxes_input, pas=pas,
                ptypes=ptypes, timeidx=timeidx, solver=solver
            )
            if solved_states == false
                solved_states = zeros(length(ele.nameinfo[:state]), length(inputs), length(timeidx))
            end
            fluxes_input = cat(fluxes_input, solved_states, dims=1)
        else
            solved_states = nothing
        end
        fluxes_outputs = run_multi_fluxes(ele, input=fluxes_input, pas=pas, ptypes=ptypes)
        fluxes_outputs_perm = permutedims(fluxes_outputs, (1, 3, 2))
        if isnothing(solved_states)
            fluxes = cat(fluxes, fluxes_outputs_perm, dims=1)
        else
            fluxes = cat(fluxes, solved_states, fluxes_outputs_perm, dims=1)
        end
    end
    if output_type == :namedtuple
        return [NamedTuple{Tuple(get_var_names(unit))}(eachslice(fluxes[:, i, :], dims=1)) for i in 1:length(inputs)]
    elseif output_type == :array
        return fluxes
    end
end
