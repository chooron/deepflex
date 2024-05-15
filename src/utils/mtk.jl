
# mtk.jl utils
function init_var_param(
    input_names::AbstractVector{Symbol},
    output_names::AbstractVector{Symbol},
    state_names::AbstractVector{Symbol},
    param_names::AbstractVector{Symbol},
)
    var_names = unique(vcat(input_names, output_names, state_names))
    varinfo = namedtuple(var_names, [first(@variables $nm(t) = 0.0) for nm in var_names])
    paraminfo = namedtuple(param_names, [first(@parameters $nm = 0.0 [tunable = true]) for nm in param_names])
    varinfo, paraminfo
end

function build_ele_system(
    funcs::AbstractVector{<:AbstractFlux},
    dfuncs::AbstractVector{<:AbstractFlux},
    varinfo::NamedTuple,
    paraminfo::NamedTuple;
    name::Symbol
)
    eqs = Equation[]
    sub_sys = ODESystem[]
    for func in funcs
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

macro itpfn(name, data, time)
    fn = Symbol(name, :_itp)
    quote
        $(esc(fn))(t) = LinearInterpolation($(esc(data)), $(esc(time)))(t)
        $(esc(fn))(t::Num) = SymbolicUtils.term($(esc(fn)), Symbolics.unwrap(t))
        $(esc(fn))(t)
    end
end

function build_itp_system(
    input::NamedTuple,
    time::AbstractVector,
    varinfo::NamedTuple;
    name::Symbol
)
    eqs = Equation[]
    for (nm, ip) in pairs(input)
        push!(eqs, varinfo[nm] ~ @itpfn(nm, ip, time))
    end
    ODESystem(eqs, t; name=Symbol(name, :itp_sys))
end

function setup_input(
    system::ODESystem;
    input::NamedTuple,
    time::AbstractVector,
    nfunc_ntp::NamedTuple,
    name::Symbol,
)
    #* 构建data的插值系统
    itp_eqs = Equation[getproperty(system, key) ~ @itpfn(key, input[key], time) for key in keys(input)]
    compose_sys = compose(ODESystem(itp_eqs, t; name=Symbol(name, :comp_sys)), system)
    sys = structural_simplify(compose_sys)
    build_u0 = Pair[]
    for k in keys(nfunc_ntp)
        func_nn_sys = getproperty(system, k)
        push!(build_u0, getproperty(getproperty(func_nn_sys, :input), :u) => zeros(length(nfunc_ntp[k])))
    end
    prob = ODEProblem(sys, build_u0, (time[1], time[end]), [], warn_initialize_determined=true)
    prob
end

function get_mtk_initstates(
    system::ODESystem;
    sol_0::ComponentVector,
    state_names::Vector{Symbol},
    nfunc_ntp::NamedTuple,
)
    #* setup init states
    u0 = [getproperty(system, nm) => sol_0[nm] for nm in keys(sol_0) if nm in state_names]
    for k in keys(nfunc_ntp)
        func_nn_sys = getproperty(system, k)
        u0 = vcat(u0, [getproperty(getproperty(func_nn_sys, :input), :u)[idx] => sol_0[nm] for (idx, nm) in enumerate(nfunc_ntp[k])])
    end
    u0
end

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