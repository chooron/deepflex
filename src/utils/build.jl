function build_flux_func(
    inputs::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
    assign_list = Assignment.(outputs, exprs)
    outputs_arr = MakeArray(outputs, Vector)
    func_args = [
        #* argument 1: Function calculation parameters
        DestructuredArgs(inputs),
        #* argument 2: Function neuralnetwork parameters
        DestructuredArgs(params),
        #* argument 3: Function time variable
        DestructuredArgs(t)
    ]
    call_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, outputs_arr, false)))
    )
    call_func
end

function build_rflux_func(
    input::Num,
    states::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
    assign_list = Assignment.(outputs, exprs)
    outputs_arr = MakeArray(outputs, Vector)
    func_args = [
        DestructuredArgs([input]),
        #* argument 1: Function calculation parameters
        DestructuredArgs([states]),
        #* argument 2: Function neuralnetwork parameters
        DestructuredArgs([params]),
        #* argument 3: Function time variable
        DestructuredArgs(t)
    ]
    call_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, outputs_arr, false)))
    )
    call_func
end

#* 构建bucket的函数
function build_ele_func(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractStateFlux},
    meta::HydroMeta
)
    #* prepare variables namedtuple
    funcs_vars = reduce(union, get_all_vars.(vcat(funcs, dfuncs)))
    funcs_vars_names = Symbolics.tosymbol.(funcs_vars, escape=false)
    funcs_vars_ntp = NamedTuple{Tuple(funcs_vars_names)}(funcs_vars)
    #* prepare variables for function building
    funcs_inputs = collect(funcs_vars_ntp[meta.inputs])
    funcs_states = collect(funcs_vars_ntp[meta.states])
    funcs_params = reduce(union, get_param_vars.(funcs))
    funcs_nns = reduce(union, get_nnparam_vars.(funcs))
    funcs_nns_bounds = nothing
    if !isempty(funcs_nns)
        funcs_nns_len = length.(funcs_nns)
        start_indices = [1; cumsum(funcs_nns_len)[1:end-1] .+ 1]
        end_indices = cumsum(funcs_nns_len)
        funcs_nns_bounds = [start:stop for (start, stop) in zip(start_indices, end_indices)]
    end

    #* Call the method in SymbolicUtils.jl to build the Flux Function
    assign_list = Assignment[]
    #* get all output flux
    output_list = Num[]

    for func in funcs
        if func isa AbstractNeuralFlux
            #* For NeuralFlux, use nn_input to match the input variable
            push!(assign_list, Assignment(func.nninfos[:inputs], MakeArray(func.inputs, Vector)))
            #* For NeuralFlux, use nn_output to match the calculation result of nn expr
            push!(assign_list, Assignment(func.nninfos[:outputs], get_exprs(func)[1]))
            #* According to the output variable name, match each index of nn_output
            for (idx, output) in enumerate(get_output_vars(func))
                push!(assign_list, Assignment(output, func.nninfos[:outputs][idx]))
                push!(output_list, output)
            end
        else
            #* According to the output variable name, match each result of the flux exprs
            for (output, expr) in zip(get_output_vars(func), get_exprs(func))
                push!(assign_list, Assignment(output, expr))
                push!(output_list, output)
            end
        end
    end
    #* convert flux output to array
    # flux_output_array = MakeArray(output_list, Vector)
    combine_list = vcat(funcs_states, output_list)
    flux_output_array = MakeArray(combine_list, Vector)

    #* Set the input argument of ODE Function
    func_args = [
        #* argument 1: Function input variables
        DestructuredArgs(funcs_inputs, :inputs, inbounds=true), 
        #* argument 2: Function state variables
        DestructuredArgs(funcs_states, :states, inbounds=true),
        #* argument 2: Function calculation parameters
        DestructuredArgs(funcs_params, :params, inbounds=true),
        #* argument 3: Function neuralnetwork parameters
        DestructuredArgs(funcs_nns, :nns, inds=funcs_nns_bounds, inbounds=true),
    ]

    #* Construct Flux Function: Func(args, kwargs, body), where body represents the matching formula between each variable and expression
    merged_flux_func = @RuntimeGeneratedFunction(toexpr(Func(func_args, [], Let(assign_list, flux_output_array, false))))

    if !isempty(dfuncs)
        #* convert diff state output to array
        diffst_output_array = MakeArray(reduce(vcat, get_exprs.(dfuncs)), Vector)
        merged_state_func = @RuntimeGeneratedFunction(
            toexpr(Func(func_args, [], Let(assign_list, diffst_output_array, false)))
        )
    else
        merged_state_func = nothing
    end
    merged_flux_func, merged_state_func
end
