using ModelingToolkit
using ModelingToolkit: t_nounits as t


macro makefn(prefix, name, args...)
    fn = Symbol(prefix * "_" * name)
    argstup = Tuple(args)
    quote
        function $(esc(fn))($(map(esc, argstup)...))
            sum([$(map(esc, argstup)...)])
        end
    end
end

macro itpfn(name, data, time)
    fn = Symbol(name, :_itp)
    name = Symbol(name)
    @variables $name(t)
    quote
        $(esc(fn))(t) = LinearInterpolation($data, $time)(t)
        @register_symbolic $(esc(fn))(t)
        $(esc(fn))(t)
    end
end

func = @itpfn("a", rand(100), 1:100)
