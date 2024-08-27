"""
$(TYPEDEF)
A struct representing common hydrological fluxes
# Fields
$(FIELDS)
# Example
```
# define a common function
function flow_func(
    i::namedtuple(:baseflow, :surfaceflow),
    p::NamedTuple;
    kw...
)
    i[:baseflow] .+ i[:surfaceflow]
end

flow_flux = SimpleFlux(
    [:baseflow, :surfaceflow],
    :flow,
    param_names=Symbol[],
    func=flow_func
)
```
"""
struct SimpleFlux <: AbstractSimpleFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "flux expressions to descripe the formula of the output variable"
    exprs::Vector{Num}
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function SimpleFlux(
        inputs::Vector{Num},
        outputs::Vector{Num},
        params::Vector{Num},
        exprs::Vector{Num},
        infos::NamedTuple
    )
        flux_func = build_flux_func(inputs, outputs, params, exprs)

        return new(
            inputs,
            outputs,
            params,
            exprs,
            flux_func,
            infos
        )
    end

    function SimpleFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        param_names::Vector{Symbol}=Symbol[];
        flux_funcs::Vector{<:Function}=Function[],
    )
        #* Get input and output names
        input_names, output_names = flux_names[1], flux_names[2]
        infos = (input=input_names, output=output_names, param=param_names, nn=Symbol[])
        if length(flux_funcs) > 0
            #* Create variables by names
            inputs = [first(@variables $var) for var in input_names]
            outputs = [first(@variables $var) for var in output_names]
            params = [first(@parameters $var) for var in param_names]
            #* When a calculation function is provided, exprs are constructed based on the calculation function and variables
            exprs = [flux_func(inputs, params) for flux_func in flux_funcs]
        else
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            inputs, outputs, params = hydro_equation.inputs, hydro_equation.outputs, hydro_equation.params
            exprs = HydroEquations.expr(hydro_equation)
        end

        #* Building the struct
        return SimpleFlux(
            inputs,
            outputs,
            params,
            exprs,
            infos
        )
    end

    function SimpleFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        exprs::Vector{Num}=Num[]
    )
        #* Get input and output variables
        inputs, outputs = fluxes[1], fluxes[2]

        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(inputs, escape=false)
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        param_names = length(params) > 0 ? Symbolics.tosymbol.(params) : Symbol[]
        infos = (input=input_names, output=output_names, param=param_names, nn=Symbol[])

        if length(exprs) == 0
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydro_equation = HydroEquation(input_names, output_names, param_names)
            exprs = HydroEquations.expr(hydro_equation)
        end

        return SimpleFlux(
            inputs,
            outputs,
            params,
            exprs,
            infos
        )
    end
end

#* callable functions
function (flux::AbstractSimpleFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
    flux.func(input, collect([pas[:params][nm] for nm in flux.infos[:param]]), nothing)
end

function (flux::AbstractSimpleFlux)(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    params_vec = collect([pas[:params][nm] for nm in flux.infos[:param]])
    reduce(hcat, flux.func.(eachslice(input, dims=2), Ref(params_vec), Ref(nothing)))
end

struct StateFlux <: AbstractStateFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of state names (Symbol) and its variables (Num)"
    state::Num
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "flux expressions to descripe the formula of the state variable"
    expr::Num
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function StateFlux(
        fluxes::Vector{Num},
        state::Num,
        params::Vector{Num}=Num[];
        expr::Num
    )
        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(fluxes, escape=false)
        state_name = Symbolics.tosymbol(state, escape=false)
        param_names = length(params) > 0 ? Symbol.(Symbolics.tosymbol.(params, escape=false)) : Symbol[]
        infos = (input=input_names, state=state_name, param=param_names, nn=Symbol[])
        state_func = build_flux_func(fluxes, [state], params, [expr])
        return new(
            fluxes,
            state,
            params,
            expr,
            state_func,
            infos
        )
    end

    function StateFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        state::Num;
    )
        influxes, outfluxes = fluxes[1], fluxes[2]
        #* Construct the default calculation formula: sum of input variables minus sum of output variables
        state_expr = sum(influxes) - sum(outfluxes)
        return StateFlux(vcat(influxes, outfluxes), state, expr=state_expr)
    end

    function StateFlux(
        states::Pair{Num,Num};
    )
        ori_state, new_state = states[1], states[2]
        #* Construct the default calculation formula: new state variable minus old state variable
        state_expr = new_state - ori_state
        return StateFlux([new_state], ori_state, expr=state_expr)
    end

    function StateFlux(
        flux_names::Vector{Symbol},
        state_name::Symbol,
        param_names::Vector{Symbol};
        state_func::Function
    )
        #* Create variables by names
        fluxes = [first(@variables $nm) for nm in flux_names]
        state = first(@variables $state_name)
        params = [first(@parameters $nm) for nm in param_names]
        state_expr = state_func(fluxes, params)
        return StateFlux(fluxes, state, params, expr=state_expr)
    end

    function StateFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        state_name::Symbol;
    )
        influx_names, outflux_names = flux_names[1], flux_names[2]
        state_input_names = vcat(influx_names, outflux_names)
        #* Construct the default calculation function: sum of input variables minus sum of output variables
        state_func = (i, p) -> sum(i[1:length(influx_names)]) - sum(i[length(influx_names)+1:length(influx_names)+length(outflux_names)])
        return StateFlux(state_input_names, state_name, Symbol[], state_func=state_func)
    end

    function StateFlux(
        state_names::Pair{Symbol,Symbol}
    )
        ori_state_name, new_state_name = state_names[1], state_names[2]
        #* Create variables by names
        ori_state = first(@variables $ori_state_name)
        new_state = first(@variables $new_state_name)
        return StateFlux(ori_state => new_state)
    end
end

struct RouteFlux{rtype} <: AbstractRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    states::Vector{Num}
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function RouteFlux(
        input::Num,
        params::Vector{Num},
        states::Vector{Num},
        routetype::Symbol,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_names = Symbolics.tosymbol.(params, escape=false)
        state_names = Symbolics.tosymbol.(states, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
        #* Setup the name information of the hydroroutement
        infos = (input=[input_name], output=[output_name], param=param_names, state=state_names)

        return new{routetype}(
            [input],
            [output],
            params,
            states,
            infos
        )
    end
end

(::AbstractRouteFlux)(::AbstractVector, ::ComponentVector) = @error "Abstract RouteFlux is not support for single timepoint"

(route::AbstractRouteFlux)(input::AbstractMatrix, pas::ComponentVector; kwargs...) = @error "Must be implemented by subtype"

function run_multi_fluxes(
    route::AbstractRouteFlux;
    input::AbstractArray,
    params::ComponentVector,
    timeidx::AbstractVector,
    kwargs...
)
    ptypes = get(kwargs, :ptypes, collect(keys(params)))
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    pytype_params = [params[ptype] for ptype in ptypes]
    sols = map(eachindex(ptypes)) do (idx)
        node_sols = reduce(hcat, route.routefunc.(eachslice(input[idx, :, :], dims=1), pytype_params[idx], Ref(timeidx)))
        node_sols
    end
    if length(sols) > 1
        sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
        return permutedims(sol_arr, (3, 1, 2))
    else
        return reshape(sols[1], 1, size(input)[3], size(input)[2])
    end
end

struct UnitHydroFlux{solvetype} <: AbstractRouteFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "unit hydrology function"
    uhfunc::Function
    "bucket information: keys contains: input, output, param"
    infos::NamedTuple

    function UnitHydroFlux(
        input::Num,
        param::Num,
        uhfunc::Function;
        solvetype::Symbol=:unithydro1,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_name = Symbolics.tosymbol(param, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
        #* Setup the name information of the hydroroutement
        infos = (input=[input_name], output=[output_name], param=[param_name], state=Symbol[])

        return new{solvetype}(
            [input],
            [output],
            [param],
            uhfunc,
            infos
        )
    end
end

"""
* 一种是用构建Discrete problem的方式求解
"""
function (flux::UnitHydroFlux{:unithydro1})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    #* 首先将lagflux转换为discrete problem
    function lag_prob(u, p, t)
        u = circshift(u, -1)
        u[end] = 0.0
        input_vec[Int(t)] .* p[:weight] .+ u
    end

    uh_weight = flux.uhfunc(pas[:params][flux.infos[:param][1]])
    prob = DiscreteProblem(lag_prob, input_vec[1] .* uh_weight, (1, length(input_vec)), ComponentVector(weight=uh_weight))
    #* 求解这个问题
    sol = solve(prob, FunctionMap())
    reshape(Array(sol)[1, :], 1, length(input_vec))
end

"""
* 一种是用构建稀疏矩阵的方式求和
"""
function (flux::UnitHydroFlux{:unithydro2})(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    uh_weight = flux.uhfunc(pas[:params][flux.infos[:param][1]])
    uh_result = [-(i - 1) => uh_wi .* input_vec for (i, uh_wi) in enumerate(uh_weight)]
    uh_sparse_matrix = spdiagm(uh_result...)
    sum_route = sum(uh_sparse_matrix, dims=2)[1:end-length(uh_weight)+1]
    reshape(sum_route, 1, length(input_vec))
end


"""
$(TYPEDEF)
A hydrological flux calculated via a neural network (based on `Lux.jl`)
# Fields
$(FIELDS)
# Example
```
et_ann = Lux.Chain(
    Lux.Dense(3 => 16, Lux.tanh),
    Lux.Dense(16 => 1, Lux.leakyrelu)
)
etnn_flux = NeuralFlux([:norm_snw, :norm_slw, :norm_temp], :evap, param_names=:etnn, chain=et_ann)
```
"""
struct NeuralFlux <: AbstractNeuralFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of output names (Symbol) and its variables (Num)"
    outputs::Vector{Num}
    "A map of neural network names (Symbol) and its variables (Num)"
    nnparam::Symbolics.Arr
    "flux expressions to descripe the formula of the output variable"
    expr::Symbolics.Arr{Num,1}
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "neural network  information: keys contains: input, output, param"
    infos::NamedTuple
    "neural network build information: keys contains: input, output, param"
    nnios::NamedTuple

    function NeuralFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        chain::Lux.AbstractExplicitContainerLayer,
    )
        #* Get input and output variables
        input_vars, output_vars = fluxes[1], fluxes[2]
        #* Get the neural network name (neural flux param name) and object
        chain_name = chain.name
        #* Initialize model parameter type for model parameter dimension definition
        init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
        init_params_axes = getaxes(init_params)

        #* Define parameter variables according to initialization parameters: Define type as Vector{parameter length}
        chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(init_params))
        #* Use Symbolics.array_term to define the slow-building parameter variables:
        #* when the model is called, it is rebuilt into the ComponentVector type according to
        #* the axes of `init_params` and the Vector type of the parameter as the calculation parameter input
        lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, init_params_axes, size=size(chain_params))

        #* Convert to a symbol based on the variable
        input_names = Symbolics.tosymbol.(input_vars, escape=false)
        output_names = Symbolics.tosymbol.(output_vars, escape=false)

        #* Constructing neural network input and output variables
        nn_input_name = Symbol(chain_name, :_input)
        nn_output_name = Symbol(chain_name, :_output)

        #* The input and output of the model can only be of type Symbolics.Arr{Num, 1},
        #* so it cannot be defined based on input_vars and output_vars
        nn_input = first(@variables $(nn_input_name)[1:length(input_names)])
        nn_output = first(@variables $(nn_output_name)[1:length(output_names)])

        #* Constructing a calculation expression based on a neural network
        flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)

        nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, init_params_axes))

        #* neuralflux infos
        infos = (input=input_names, output=output_names, param=Symbol[], nn=[chain_name])
        new(
            input_vars,
            output_vars,
            chain_params,
            flux_expr,
            nn_func,
            infos,
            (input=nn_input, output=nn_output, paramlen=length(init_params)),
        )
    end

    function NeuralFlux(
        flux_names::Pair{Vector{Symbol},Vector{Symbol}},
        chain::Lux.AbstractExplicitContainerLayer,
    )
        input_names, output_names = flux_names[1], flux_names[2]

        input_vars = [first(@variables $input_name) for input_name in input_names]
        output_vars = [first(@variables $output_name) for output_name in output_names]

        return NeuralFlux(input_vars => output_vars, chain)
    end
end

function (flux::AbstractNeuralFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][flux.infos[:nn][1]]
    flux.func(input, nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::AbstractMatrix, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][flux.infos[:nn][1]]
    flux.func(input', nn_params_vec)
end
