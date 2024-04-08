# mtk.jl utils
function init_var_param(nameinfo::NameInfo)
    var_names = vcat(nameinfo.input_names, nameinfo.output_names, nameinfo.state_names)
    varinfo = namedtuple(var_names, [first(@variables $nm(t)) for nm in var_names])
    paraminfo = namedtuple(nameinfo.param_names, [first(@parameters $nm) for nm in nameinfo.param_names])
    varinfo, paraminfo
end

function build_ele_system(
    funcs::AbstractVector{AbstractFlux},
    dfuncs::AbstractVector{AbstractFlux},
    varinfo::NamedTuple,
    paraminfo::NamedTuple;
    name::Symbol
)
    eqs = Equation[]
    for func in funcs
        tmp_input = namedtuple(get_input_names(func), [varinfo[nm] for nm in get_input_names(func)])
        tmp_param = namedtuple(get_param_names(func), [paraminfo[nm] for nm in get_param_names(func)])
        push!(eqs, varinfo[first(get_output_names(func))] ~ func(tmp_input, tmp_param)[first(get_output_names(func))])
    end
    for dfunc in dfuncs
        tmp_input = namedtuple(get_input_names(dfunc), [varinfo[nm] for nm in get_input_names(dfunc)])
        tmp_param = namedtuple(get_param_names(dfunc), [paraminfo[nm] for nm in get_param_names(dfunc)])
        push!(eqs, D(varinfo[first(get_output_names(dfunc))]) ~ dfunc(tmp_input, tmp_param)[first(get_output_names(dfunc))])
    end
    ODESystem(eqs, t; name=Symbol(name, :sys))
end

function build_itp_system(
    input::NamedTuple,
    time::AbstractVector,
    varinfo::NamedTuple;
    name::Symbol
)
    eqs = Equation[]
    for nm in keys(input)
        tmp_itp = LinearInterpolation(input[nm], time, extrapolate=true)
        push!(eqs, varinfo[nm] ~ tmp_itp(t))
    end
    ODESystem(eqs, t; name=Symbol(name, :itp_sys))
end