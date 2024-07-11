function build_state_func(
    funcs::Vector{<:AbstractFlux},
    dfunc::AbstractStateFlux,
    input_names::Vector{Symbol},
)
    fluxes_vars_ntp = reduce(merge, [merge(func.input_info, func.output_info) for func in vcat(funcs, [dfunc])])
    funcs_params_ntp = reduce(merge, [func.param_info for func in vcat(funcs, [dfunc])])
    #* 构建state计算的函数并将所有中间状态替换
    substitute_vars_dict = Dict()
    for var_nm in keys(fluxes_vars_ntp)
        for func in funcs
            tmp_output_names = get_output_names(func)
            for j in eachindex(tmp_output_names)
                if var_nm == tmp_output_names[j]
                    substitute_vars_dict[fluxes_vars_ntp[var_nm]] = func.flux_exprs[j]
                end
            end
        end
    end
    state_expr_sub = dfunc.state_expr
    for _ in 1:length(substitute_vars_dict)
        state_expr_sub = substitute(state_expr_sub, substitute_vars_dict)
    end
    state_func = build_function(state_expr_sub, collect(fluxes_vars_ntp[input_names]), collect(funcs_params_ntp), expression=Val{false})
    state_func
end

function build_state_funcv2(
    funcs::Vector{<:AbstractFlux},
    dfunc::AbstractStateFlux,
    input_names::Vector{Symbol},
)
    fluxes_vars_ntp = reduce(merge, [merge(func.input_info, func.output_info) for func in vcat(funcs, [dfunc])])
    funcs_params_ntp = reduce(merge, [func.param_info for func in vcat(funcs, [dfunc])])
    
    assign_list = Assignment[]
    for func in funcs
        if func isa AbstractNeuralFlux
            push!(assign_list, Assignment(func.nn_info[:input], MakeArray(collect(func.input_info), Vector)))
            push!(assign_list, Assignment(func.nn_info[:output], func.flux_expr))
            for (idx, output) in enumerate(collect(func.output_info))
                push!(assign_list, Assignment(output, func.nn_info[:output][idx]))
            end
        else
            for (output, expr) in zip(collect(func.output_info), func.flux_exprs)
                push!(assign_list, Assignment(output, expr))
            end
        end
    end

    func_args = [DestructuredArgs(collect(fluxes_vars_ntp[input_names])), DestructuredArgs(collect(funcs_params_ntp))]
    merged_state_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, dfunc.state_expr, false)))
    )
    merged_state_func
end

function build_state_funcv3(
    funcs::Vector{<:AbstractFlux},
    dfunc::AbstractStateFlux,
    input_names::Vector{Symbol},
)
    total_param_names = keys(reduce(hcat, [func.param_info for func in vcat(funcs, dfunc)]))
    state_input_names = keys(dfunc.input_info)
    state_param_names = keys(dfunc.param_info)

    func_input_names = [get_input_names(func) for func in funcs]
    func_output_names = [get_output_names(func) for func in funcs]
    func_param_names = [get_param_names(func) for func in funcs]

    flux_exprs = [
        :(($(output_name...),) = $flux([$(input_name...),], [$(param_name...)]))
        for (input_name, output_name, param_name, flux) in zip(func_input_names, func_output_names, func_param_names, funcs)
    ]

    state_func = build_function(dfunc.state_expr, collect(dfunc.input_info), collect(dfunc.param_info), expression=Val{true})

    state_func_expr = :((i, p) -> begin
        ($(input_names...),) = i
        ($(total_param_names...),) = p
        $(flux_exprs...)
        return $(state_func)([$(state_input_names...)], [$(state_param_names...)])
    end)

    merged_state_func = @RuntimeGeneratedFunction(state_func_expr)
    merged_state_func
end