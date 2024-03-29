using ModelingToolkit

@variables t
D = Differential(t)

function decay(; name)
    @parameters a
    @variables x(t) f(t)
    ODESystem([
            D(x) ~ -a * x + f
        ], t;
        name = name)
end

decay1 = decay(name=:decay1)
decay2 = decay(name=:decay1)

connected = compose(
    ODESystem([decay2.f ~ decay1.x,
               D(decay1.f) ~ 0], t; name = :connected), decay1, decay2)

equations(connected)

#4-element Vector{Equation}:
# Differential(t)(decay1₊f(t)) ~ 0
# decay2₊f(t) ~ decay1₊x(t)
# Differential(t)(decay1₊x(t)) ~ decay1₊f(t) - (decay1₊a*(decay1₊x(t)))
# Differential(t)(decay2₊x(t)) ~ decay2₊f(t) - (decay2₊a*(decay2₊x(t)))

simplified_sys = structural_simplify(connected)

equations(simplified_sys)