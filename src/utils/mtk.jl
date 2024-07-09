"""
$(SIGNATURES)

Based on all ordinary hydrological calculation fluxes and state hydrological fluxes in the hydrological calculation module,
constructing a system of ordinary differential equations based on ModelingToolkit.jl.
"""
function build_ele_system(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractFlux};
    name::Symbol,
)
    funcs_eqs = reduce(vcat, [func.flux_eqs for func in funcs])
    dfuncs_eqs = reduce(vcat, [dfunc.state_eq for dfunc in dfuncs])
    sub_sys = reduce(vcat, [nfunc.sub_sys for nfunc in filter(f -> f isa AbstractNeuralFlux, funcs)])
    base_sys = ODESystem(vcat(funcs_eqs, dfuncs_eqs), t; name=Symbol(name, :_base_sys), systems=sub_sys)
    base_sys
end

function build_state_func(
    fluxes::Vector{<:AbstractFlux},
    state_expr::Num,
    fluxes_vars_ntp::NamedTuple,
    funcs_params::NamedTuple,
    state_input_names::Vector{Symbol},
)
    #* 构建state计算的函数并将所有中间状态替换
    substitute_vars_dict = Dict()
    for var_nm in keys(fluxes_vars_ntp)
        for flux in fluxes
            tmp_output_names = get_output_names(flux)
            for j in eachindex(tmp_output_names)
                if var_nm == tmp_output_names[j]
                    tmp_flux_exprs = Symbolics.rhss(flux.flux_eqs)
                    substitute_vars_dict[fluxes_vars_ntp[var_nm]] = tmp_flux_exprs[j]
                end
            end
        end
    end
    state_expr_sub = state_expr
    for _ in 1:length(substitute_vars_dict)
        state_expr_sub = substitute(state_expr_sub, substitute_vars_dict)
    end
    state_func = build_function(state_expr_sub, collect(fluxes_vars_ntp[state_input_names]), collect(funcs_params), expression=Val{false})
    state_func
end

function build_state_funcv2(
    fluxes::Vector{<:AbstractFlux},
    state_func::Function,
    state_input_names::Vector{Symbol},
    state_param_names::Vector{Symbol},
)
    flux_exprs = [
        quote
            ($(get_output_names(flux)...),) = $(flux)([$(get_input_names(flux))...], [$(get_param_names(flux))...])
        end
        for flux in fluxes
    ]

    total_param_names = union(state_param_names, get_param_names(fluxes))

    state_func_expr = :((i, p) -> begin
        ($(state_input_names...),) = i
        ($(total_param_names...),) = p
        $flux_exprs
        return $(state_func)([$(state_input_names)...], [$(state_param_names)...])
    end)

    func = @RuntimeGeneratedFunction(state_func_expr)
    func
end

"""
$(SIGNATURES)

A macro for building and registering interpolation functions
"""
macro itpfn(name, data, ts)
    fn = Symbol(name, :_itp)
    quote
        $(esc(fn))(t) = LinearInterpolation($(esc(data)), $(esc(ts)))(t)
        $(esc(fn))(t::Num) = SymbolicUtils.term($(esc(fn)), Symbolics.unwrap(t))
        $(esc(fn))(t)
    end
end

"""
$(SIGNATURES)

Build an interpolation system through input data and couple it with the system 
initialized by hydrological element to form a complete system
"""
function setup_input(
    ele_system::ODESystem,
    input::NamedTuple,
    timeidx::Vector,
    input_names::Vector,
    name::Symbol,
)
    #* 构建data的插值系统
    #! 这个语句会影响zygote的优化计算
    itp_eqs = Equation[getproperty(ele_system, nm) ~ @itpfn(nm, input[nm], timeidx) for nm in input_names]
    compose_sys = complete(ODESystem(itp_eqs, t; name=Symbol(name, :comp_sys), systems=[ele_system]))
    structural_simplify(compose_sys)
end

"""
$(SIGNATURES)

Used to construct the initialization problem based on the complete system, after constructing the complete system,
taking into account that systems with neurohydrological fluxes often have special requirements for problem initialization
used for hydrological element
"""
function init_prob(
    build_system::ODESystem,
    prebuild_system::ODESystem,
    nfunc_ntp::NamedTuple,
    time::AbstractVector,
)
    #* 设置系统初始值并构建problem
    build_u0 = Pair[]
    for k in keys(nfunc_ntp)
        func_nn_sys = getproperty(prebuild_system, k)
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(nfunc_ntp[k])))
    end
    prob = ODEProblem(build_system, build_u0, (time[1], time[end]), Pair[], warn_initialize_determined=false)
    prob
end

"""
$(SIGNATURES)

Obtain initialization state for problem `remake`
based on the re-entered initial state combined with the constructed problem
"""
function get_mtk_initstates(
    prebuild_systems::Vector{ODESystem},
    sol_0::NamedTuple,
    ele_state_names::Vector,
    nfunc_ntps::Vector,
)
    #* setup init states
    u0 = Pair[]
    for (prebuild_system, state_names, nfunc_ntp) in zip(prebuild_systems, ele_state_names, nfunc_ntps)
        for nm in state_names
            push!(u0, getproperty(prebuild_system, nm) => sol_0[nm])
        end
        for k in keys(nfunc_ntp)
            func_nn_sys = getproperty(prebuild_system, k)
            for (idx, nm) in enumerate(nfunc_ntp[k])
                push!(u0, getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm])
            end
        end
    end
    u0
end

"""
$(SIGNATURES)

Obtain update parameters for problem `remake`
based on the re-entered parameters combined with the constructed problem
"""
function get_mtk_params(
    prebuild_systems::Vector{ODESystem},
    params::NamedTuple,
)
    #* setup parameters
    p = Pair[]
    for prebuild_system in prebuild_systems
        for nm in ModelingToolkit.parameters(prebuild_system)
            if contains(string(nm), "₊")
                tmp_nn = split(string(nm), "₊")[1]
                push!(p, getproperty(getproperty(prebuild_system, Symbol(tmp_nn)), :p) => Vector(params[Symbol(tmp_nn)]))
            else
                push!(p, getproperty(prebuild_system, Symbol(nm)) => params[Symbol(nm)])
            end
        end
    end
    p
end
