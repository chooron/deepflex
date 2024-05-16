"""
$(SIGNATURES)

Init system variables based on provided variable names (including input, output, state names),
and Init system parameters based on parameter names, which are used in the initialization 
of the hydrological calculation module
"""
function init_var_param(
    var_names::Vector{Symbol},
    param_names::Vector{Symbol},
)
    varinfo = namedtuple(var_names, [first(@variables $nm(t) = 0.0) for nm in var_names])
    paraminfo = namedtuple(param_names, [first(@parameters $nm = 0.0 [tunable = true]) for nm in param_names])
    varinfo, paraminfo
end

"""
$(SIGNATURES)

Based on all ordinary hydrological calculation fluxes and state hydrological fluxes in the hydrological calculation module,
constructing a system of ordinary differential equations based on ModelingToolkit.jl.
"""
function build_ele_system(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractFlux};
    var_names::Vector{Symbol},
    param_names::Vector{Symbol},
    name::Symbol
)
    varinfo, paraminfo = init_var_param(var_names, param_names)
    eqs = Equation[]
    sub_sys = ODESystem[]
    for func in funcs
        #* Here is a system build for adaptation to neurohydrological fluxes
        if func isa NeuralFlux
            nn_in = RealInputArray(nin=length(func.input_names), name=Symbol(func.chain_name, :_in_sys))
            nn_out = RealOutputArray(nout=length(func.output_names), name=Symbol(func.chain_name, :_out_sys))
            eqs = vcat(eqs, [varinfo[nm] ~ nn_in.u[idx] for (idx, nm) in enumerate(get_input_names(func))])
            eqs = vcat(eqs, [varinfo[nm] ~ nn_out.u[idx] for (idx, nm) in enumerate(get_output_names(func))])
            eqs = vcat(eqs, [connect(nn_in, func.nn.input), connect(nn_out, func.nn.output)])
            sub_sys = vcat(sub_sys, [func.nn, nn_in, nn_out])
        else
            tmp_input = namedtuple(get_input_names(func), [varinfo[nm] for nm in get_input_names(func)])
            tmp_param = namedtuple(get_param_names(func), [paraminfo[nm] for nm in get_param_names(func)])
            eqs = vcat(eqs, [varinfo[nm] ~ func(tmp_input, tmp_param)[nm] for nm in get_output_names(func)])
        end
    end
    for dfunc in dfuncs
        tmp_input = namedtuple(get_input_names(dfunc), [varinfo[nm] for nm in get_input_names(dfunc)])
        tmp_param = namedtuple(get_param_names(dfunc), [paraminfo[nm] for nm in get_param_names(dfunc)])
        eqs = vcat(eqs, [D(varinfo[nm]) ~ dfunc(tmp_input, tmp_param)[nm] for nm in get_output_names(dfunc)])
    end
    base_sys = ODESystem(eqs, t; name=Symbol(name, :base_sys), systems=sub_sys)
    base_sys
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
    system::ODESystem;
    input::NamedTuple,
    time::Vector,
    name::Symbol,
)
    #* 构建data的插值系统
    itp_eqs = Equation[getproperty(system, key) ~ @itpfn(key, input[key], time) for key in keys(input)]
    compose_sys = compose(ODESystem(itp_eqs, t; name=Symbol(name, :comp_sys)), system)
    structural_simplify(compose_sys)
end

"""
$(SIGNATURES)

Build an interpolation system through input data and couple it with the system 
initialized by hydrological unit to form a complete system
"""
function setup_input(
    unit_system::ODESystem,
    ele_systems::Vector{ODESystem};
    input::NamedTuple,
    input_names::Vector{Vector{Symbol}},
    name::Symbol,
)
    eqs = Equation[]
    for (idx, system) in enumerate(ele_systems)
        for nm in filter(nm -> nm in input_names[idx], keys(input))
            push!(eqs, getproperty(system, nm) ~ @itpfn(nm, input[nm], input[:time]))
        end
    end
    compose_sys = compose(ODESystem(eqs, t; name=Symbol(name, :_comp_sys)), unit_system)
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
    system::ODESystem;
    nfunc_ntp::NamedTuple,
    time::Vector,
)
    #* 设置系统初始值并构建problem
    build_u0 = Pair[]
    for k in keys(nfunc_ntp)
        func_nn_sys = getproperty(system, k)
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(nfunc_ntp[k])))
    end
    prob = ODEProblem(build_system, build_u0, (time[1], time[end]), [], warn_initialize_determined=true)
    prob
end

"""
$(SIGNATURES)

Used to construct the initialization problem based on the complete system, after constructing the complete system,
taking into account that systems with neurohydrological fluxes often have special requirements for problem initialization
used for hydrological unit
"""
function init_prob(
    build_system::ODESystem,
    systems::Vector{ODESystem};
    nfunc_ntp_list::Vector,
    time::Vector
)
    #* 设置系统初始值并构建problem
    build_u0 = Pair[]
    for (idx, system) in enumerate(systems)
        for k in keys(nfunc_ntp_list[idx])
            func_nn_sys = getproperty(system, k)
            push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(nfunc_ntp_list[idx][k])))
        end
    end
    prob = ODEProblem(build_system, build_u0, (time[1], time[end]), [], warn_initialize_determined=true)
    prob
end

"""
$(SIGNATURES)

Obtain initialization state for problem `remake`
based on the re-entered initial state combined with the constructed problem
"""
function get_mtk_initstates(
    system::ODESystem;
    sol_0::NamedTuple,
    state_names::Vector{Symbol},
    nfunc_ntp::NamedTuple,
)
    #* setup init states
    u0 = Pair[getproperty(system, nm) => sol_0[nm] for nm in keys(sol_0) if nm in state_names]
    for k in keys(nfunc_ntp)
        func_nn_sys = getproperty(system, k)
        u0 = vcat(u0, Pair[getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm] for (idx, nm) in enumerate(nfunc_ntp[k])])
    end
    u0
end
"""
$(SIGNATURES)

Obtain update parameters for problem `remake`
based on the re-entered parameters combined with the constructed problem
"""
function get_mtk_params(
    system::ODESystem;
    params::ComponentVector,
)
    #* setup init states
    p = Pair[]
    for nm in ModelingToolkit.parameters(system)
        if contains(string(nm), "₊")
            tmp_nn = split(string(nm), "₊")[1]
            push!(p, getproperty(getproperty(system, Symbol(tmp_nn)), :p) => Vector(params[Symbol(tmp_nn)]))
        else
            push!(p, getproperty(system, Symbol(nm)) => params[Symbol(nm)])
        end
    end
    p
end