"""
    HydroModel <: AbstractModel

Represents a hydrological model composed of multiple components.

# Fields
- `infos::NamedTuple`: Contains metadata about the model, including name, input variables, all variables, output variables, state variables, and neural network variables.
- `components::Vector{<:AbstractComponent}`: A vector of hydrological computation elements (components) that make up the model.
- `varindices::Vector`: A vector of indices for each component's input, used to map overall model inputs to component-specific inputs.

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

Each component's kwargs may be different, include solver, interp
"""
struct HydroModel <: AbstractModel
    "meta data of hydrological model"
    meta::HydroMeta
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "input variables index for each components"
    varindices::Vector
    "all variables names"
    var_names::Vector{Symbol}

    function HydroModel(;name::Symbol, components::Vector{<:AbstractComponent})
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
        model_meta = HydroMeta(
            name=name,
            inputs=input_names, outputs=output_names,
            states=state_names, params=param_names,
            nns=nn_names
        )
        new(
            model_meta,
            components,
            input_idx,
            var_names,
        )
    end
end

function (model::HydroModel)(
    input::NamedTuple,
    pas::ComponentVector,
    timeidx::Vector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)   
    @assert all(nm -> nm in keys(input), get_input_names(model)) "input must contain all input names"
    input_matrix = Matrix(reduce(hcat, [input[nm] for nm in get_input_names(model)])')
    return model(input_matrix, pas, timeidx; config=config, kwargs...)
end

# 求解并计算
function (model::HydroModel)(
    input::Matrix,
    pas::ComponentVector,
    timeidx::Vector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)
    #* 如果compkwargs是NamedTuple,则将其填充为Vector{NamedTuple}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    @assert size(input, 1) == length(get_input_names(model)) "input matrix must have the same number of columns as the input names"
    @assert size(input, 2) == length(timeidx) "input matrix must have the same number of rows as the timeidx length"
    fluxes = input
    for (comp_, idx, config_) in zip(model.components, model.varindices, comp_configs)
        # Ensure convert_to_ntp is false for inner component computation
        if comp_ isa AbstractEstimator
            tmp_pas = comp_(fluxes[idx, :], pas, timeidx; config=config_, convert_to_ntp=false)
            pas = update_ca(pas, tmp_pas)
        else
            tmp_fluxes = comp_(fluxes[idx, :], pas, timeidx; config=config_, convert_to_ntp=false)
            fluxes = cat(fluxes, tmp_fluxes, dims=1)
        end
    end
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)
    if convert_to_ntp
        return NamedTuple{Tuple(model.var_names)}(eachrow(fluxes))
    else
        return fluxes
    end
end

#* 多输入构建大型方程求解并计算
function (model::HydroModel)(
    inputs::Vector{<:NamedTuple},
    pas::ComponentVector,
    timeidx::Vector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=(solver=ODESolver(), ptypes=keys(pas[:params]), interp=LinearInterpolation),
    kwargs...
)
    #* 如果compkwargs是NamedTuple,则将其填充为Vector{NamedTuple}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"

    fluxes = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, [input[nm] for nm in get_input_names(model)]) for input in inputs])
    fluxes = permutedims(fluxes, (2, 3, 1))
    for (comp_, idx_, config_) in zip(model.components, model.varindices, comp_configs)
        if comp_ isa AbstractEstimator
            tmp_pas = comp_(fluxes[idx_, :, :], pas, timeidx; config=config_, convert_to_ntp=false)
            pas = update_ca(pas, tmp_pas)
        else
            tmp_fluxes = comp_(fluxes[idx_, :, :], pas, timeidx; config=config_, convert_to_ntp=false)
            fluxes = cat(fluxes, tmp_fluxes, dims=1)
        end
    end
    convert_to_ntp = get(kwargs, :convert_to_ntp, false)
    if convert_to_ntp
        return [NamedTuple{Tuple(model.var_names)}(eachslice(fluxes[:, i, :], dims=1)) for i in 1:length(inputs)]
    else
        return fluxes
    end
end
