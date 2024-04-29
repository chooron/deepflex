
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
            # for (idx, nm) in enumerate(get_input_names(func))
            #     push!(eqs, varinfo[nm] ~ func.nn.input.u[idx])
            # end
            # for (idx, nm) in enumerate(get_output_names(func))
            #     push!(eqs, varinfo[nm] ~ func.nn.output.u[idx])
            # end
            # push!(sub_sys, func.nn)
            nn_in = RealInputArray(nin=length(func.input_names), name=Symbol(func.param_names, :_in_sys))
            nn_out = RealOutputArray(nout=length(func.output_names), name=Symbol(func.param_names, :_out_sys))
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

# function build_nn_system(
#     input_dim::Int=3,
#     output_dim::Int=1,
#     chain,
#     rng,
#     eltype=Float64
# )
#     lux_p = Lux.initialparameters(rng, chain)
#     ca = ComponentArray{eltype}(lux_p)

#     @parameters p[1:length(ca)] = Vector(ca)
#     @parameters T::typeof(typeof(p)) = typeof(p) [tunable = false]

#     @named input = RealInput(nin=input_dim)
#     @named output = RealOutput(nout=output_dim)

#     out = stateless_apply(chain, input.u, lazyconvert(typeof(ca), p))

#     eqs = [output.u ~ out]

#     @named ude_comp = ODESystem(
#         eqs, t_nounits, [], [p, T], systems=[input, output])
#     return ude_comp
# end