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
    ]
    call_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, outputs_arr, false)))
    )
    call_func
end

function build_flux_func_with_time(
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

#* 构建bucket的函数
function build_ele_func(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractStateFlux},
    infos::NamedTuple
)
    #* prepare variables namedtuple
    funcs_vars = reduce(union, get_all_vars.(vcat(funcs, dfuncs)))
    funcs_vars_names = Symbolics.tosymbol.(funcs_vars, escape=false)
    funcs_vars_ntp = NamedTuple{Tuple(funcs_vars_names)}(funcs_vars)
    #* prepare variables for function building
    funcs_inputs = collect(funcs_vars_ntp[infos.input])
    funcs_states = collect(funcs_vars_ntp[infos.state])
    funcs_params = reduce(union, get_param_vars.(funcs))
    funcs_nns = reduce(union, get_nnparam_vars.(funcs))

    #* Call the method in SymbolicUtils.jl to build the Flux Function
    assign_list = Assignment[]
    #* get all output flux
    output_list = Num[]

    for func in funcs
        if func isa AbstractNeuralFlux
            #* For NeuralFlux, use nn_input to match the input variable
            push!(assign_list, Assignment(func.nnios[:input], MakeArray(func.inputs, Vector)))
            #* For NeuralFlux, use nn_output to match the calculation result of nn expr
            push!(assign_list, Assignment(func.nnios[:output], get_exprs(func)[1]))
            #* According to the output variable name, match each index of nn_output
            for (idx, output) in enumerate(get_output_vars(func))
                push!(assign_list, Assignment(output, func.nnios[:output][idx]))
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
    flux_output_array = MakeArray(output_list, Vector)

    #* Set the input argument of ODE Function
    func_args = [
        #* argument 1: Function input and state variables
        DestructuredArgs(vcat(funcs_inputs, funcs_states)),
        #* argument 2: Function calculation parameters
        DestructuredArgs(funcs_params),
        #* argument 3: Function neuralnetwork parameters
        DestructuredArgs(funcs_nns),
        #* argument 4: current time idx for time-varying flux
        DestructuredArgs(t),
    ]

    #* Construct Flux Function: Func(args, kwargs, body), where body represents the matching formula between each variable and expression
    merged_flux_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, flux_output_array, false)))
    )

    if length(dfuncs) > 0
        #* convert diff state output to array
        diffst_output_array = MakeArray(reduce(vcat, get_exprs.(dfuncs)), Vector)
        #* Set the input argument of ODE Function
        dfunc_args = [
            #* argument 1: Function input variables
            DestructuredArgs(funcs_inputs),
            #* argument 2: Function state variables
            DestructuredArgs(funcs_states),
            #* argument 3: Function calculation parameters
            DestructuredArgs(funcs_params),
            #* argument 4: Function neuralnetwork parameters
            DestructuredArgs(funcs_nns),
            #* argument 5: current time idx for time-varying flux
            DestructuredArgs(t),
        ]
        merged_state_func = @RuntimeGeneratedFunction(
            toexpr(Func(dfunc_args, [], Let(assign_list, diffst_output_array, false)))
        )
    else
        merged_state_func = nothing
    end
    merged_flux_func, merged_state_func
end
