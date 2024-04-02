
function solve_prob(
    ele::HydroElement;
    input::NamedTuple,
    params::Union{NamedTuple,ComponentVector},
    init_states::Union{NamedTuple,ComponentVector}
)
    # fit interpolation functions
    ele_input_names = ele.nameinfo.input_names
    ele_state_names = ele.nameinfo.state_names
    itp_dict = namedtuple(ele_input_names, [LinearInterpolation(input[nm], input[:time], extrapolate=true) for nm in ele_input_names])
    function singel_ele_ode_func!(du, u, p, t)
        tmp_input = namedtuple(ele_input_names, [itp_dict[nm](t) for nm in ele_input_names])
        tmp_states = namedtuple(ele_state_names, u)
        tmp_fluxes = merge(tmp_input, tmp_states)
        for func in ele.funcs
            tmp_fluxes = merge(tmp_fluxes, func(tmp_fluxes, params))
        end 
        du[:] = [dfunc(tmp_fluxes, params)[dfunc.output_names] for dfunc in ele.dfuncs]
    end
    prob = ODEProblem(singel_ele_ode_func!, collect(init_states[ele_state_names]), (input[:time][1], input[:time][end]))
    sol = solve(prob, Tsit5(), saveat=input[:time])
    solved_u = hcat(sol.u...)
    sol
    # state_names = collect(keys(init_states))
    # namedtuple(state_names, [solved_u[idx, :] for idx in 1:length(state_names)])
end