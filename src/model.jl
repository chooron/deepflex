"""
    HydroModel <: AbstractModel

Represents a hydrological model composed of multiple components.

# Fields
- `infos::NamedTuple`: Contains metadata about the model, including name, input variables, all variables, output variables, state variables, and neural network variables.
- `components::Vector{<:AbstractComponent}`: A vector of hydrological computation elements (components) that make up the model.
- `input_idx::Vector`: A vector of indices for each component's input, used to map overall model inputs to component-specific inputs.

# Constructor
    HydroModel(name; components::Vector{<:AbstractComponent})

Constructs a HydroModel with the given name and components.

# Description
HydroModel is a structure that encapsulates a complete hydrological model. It manages multiple hydrological components, 
handles the flow of data between these components, and provides methods for running simulations.

The model automatically determines the connections between components based on their input and output variables. 
It also keeps track of all variables in the system, including inputs, outputs, states, and any neural network parameters.

When called as a function, the HydroModel applies its components in sequence, passing the outputs of earlier components 
as inputs to later ones, effectively simulating the hydrological system over time.

"""
struct HydroModel <: AbstractModel
    "hydrological computation model information"
    infos::NamedTuple
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "input idx for each components"
    input_idx::Vector

    function HydroModel(name; components::Vector{<:AbstractComponent})
        #* 获取每个element的输出结果,然后与输入结果逐次拼接,获取每次输入的matrix的idx
        input_names, output_names, state_names = get_var_names(components)
        nn_names = reduce(union, get_nn_names.(components))
        param_names = reduce(union, get_param_names.(components))
        var_names = input_names
        input_idx = Vector[]
        for component in components
            tmp_input_idx = map(get_input_names(component)) do nm
                findfirst(varnm -> varnm == nm, var_names)
            end
            #* 更新model_var_names
            var_names = reduce(vcat, [var_names, get_state_names(component), get_output_names(component)])
            push!(input_idx, tmp_input_idx)
        end
        model_infos = (name=name, input=input_names, var=var_names, output=output_names, state=state_names, nn=nn_names, param=param_names)
        new(
            model_infos,
            components,
            input_idx,
        )
    end
end

# 求解并计算
function (model::HydroModel)(
    input::NamedTuple,
    pas::ComponentVector;
    timeidx::Vector,
    solver::AbstractSolver=ODESolver(),
    output_type::Symbol=:namedtuple
)
    fluxes = Matrix(reduce(hcat, [input[nm] for nm in get_input_names(model)])')
    for (tmp_comp, idx) in zip(model.components, model.input_idx)
        # todo 添加对于estimator的识别
        if tmp_comp isa AbstractEstimator
            tmp_pas = tmp_comp(fluxes[idx, :], pas, timeidx=timeidx, solver=solver)
            pas = update_ca(pas, tmp_pas)
        else
            tmp_fluxes = tmp_comp(fluxes[idx, :], pas, timeidx=timeidx, solver=solver)
            fluxes = cat(fluxes, tmp_fluxes, dims=1)
        end
    end
    if output_type == :namedtuple
        return NamedTuple{Tuple(get_var_names(model))}(eachrow(fluxes))
    elseif output_type == :matrix
        return fluxes
    else
        error("output_type must be :namedtuple or :matrix")
    end
end

#* 多输入构建大型方程求解并计算
function (model::HydroModel)(
    inputs::Vector{<:NamedTuple},
    pas::ComponentVector;
    timeidx::Vector,
    ptypes::Vector{Symbol}=collect(keys(pas[:params])),
    solver::AbstractSolver=ODESolver(),
    output_type::Symbol=:array
)
    fluxes = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, [input[nm] for nm in get_input_names(model)]) for input in inputs])
    fluxes = permutedims(fluxes, (2, 3, 1))
    for (tmp_comp, idx) in zip(model.components, model.input_idx)
        fluxes_outputs = tmp_comp(fluxes[idx, :, :], pas, ptypes=ptypes, timeidx=timeidx, solver=solver)
        fluxes = cat(fluxes, fluxes_outputs, dims=1)
    end
    if output_type == :namedtuple
        return [NamedTuple{Tuple(get_var_names(model))}(eachslice(fluxes[:, i, :], dims=1)) for i in 1:length(inputs)]
    elseif output_type == :array
        return fluxes
    else
        error("output_type must be :namedtuple or :array")
    end
end
