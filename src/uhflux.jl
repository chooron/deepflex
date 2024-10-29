struct UHFunction{uhtype}
    function UHFunction(uhtype::Symbol)
        return new{uhtype}()
    end
end

function (uh::UHFunction{:UH_1_HALF})(t, lag)
    if t - lag > 0
        typeof(lag)(1)
    else
        (t / lag)^2.5
    end
end

get_uh_tmax(::UHFunction{:UH_1_HALF}, lag) = lag

function (uh::UHFunction{:UH_2_FULL})(t, lag)
    if t - lag * 2 > 0
        typeof(lag)(1)
    elseif t - lag > 0
        (1 - 0.5 * abs(2 - t / lag)^2.5)
    else
        (0.5 * abs(t / lag)^2.5)
    end
end

get_uh_tmax(::UHFunction{:UH_2_FULL}, lag) = 2 * lag

"""
    UnitHydroFlux{solvetype} <: AbstractRouteFlux

Represents a unit hydrograph flux model for routing water through a hydrological system.

# Fields
- `inputs::Vector{Num}`: A vector of input variables (Num).
- `outputs::Vector{Num}`: A vector of output variables (Num).
- `params::Vector{Num}`: A vector of parameter variables (Num).
- `uhfunc::Function`: The unit hydrograph function.
- `meta::HydroMeta`: A named tuple containing information about inputs, outputs, parameters, and states.

# Constructor
    UnitHydroFlux(input::Num, param::Num, uhfunc::Function; solvetype::Symbol=:unithydro1)

# Arguments
- `input::Num`: The input variable.
- `param::Num`: The parameter variable.
- `uhfunc::Function`: The unit hydrograph function.
- `solvetype::Symbol`: The solver type (default is `:unithydro1`).

# Description
UnitHydroFlux represents a unit hydrograph flux model for routing water through a hydrological system. 
It uses a unit hydrograph function to transform input flows into routed output flows.

The structure supports different solving methods, specified by the `solvetype` parameter. 
Currently, it implements two solver types:
- `:unithydro1`: Uses a discrete problem approach to calculate the routed flow.
- `:unithydro2`: Uses a sparse matrix approach for more efficient computation, especially for longer time series.
The choice of solver type can affect both the performance and memory usage of the model.

This flux model is particularly useful in hydrological modeling for representing the 
temporal distribution of runoff as it moves through a watershed. It can account for the 
lag and attenuation of flow as water travels through the system.

The `uhfunc` field holds the unit hydrograph function, which typically takes a parameter 
(like time) and returns weights that describe how an input is distributed over time in the output.

When called, the UnitHydroFlux object applies the unit hydrograph to the input flow series, 
effectively convolving the input with the unit hydrograph to produce the routed output flow.

This structure is designed to be flexible and can be integrated into larger hydrological models 
to represent various routing processes in different parts of a water system.

"""
struct UnitHydroFlux{solvetype} <: AbstractUnitHydroFlux
    "A vector of input variables (Num)"
    inputs::Vector{Num}
    "A vector of output variables (Num)"
    outputs::Vector{Num}
    "A vector of parameter variables (Num)"
    params::Vector{Num}
    "The unit hydrograph function"
    uhfunc::UHFunction
    "A named tuple containing information about inputs, outputs, parameters, and states"
    meta::HydroMeta

    function UnitHydroFlux(
        input::Num,
        param::Num,
        uhtype::Symbol;
        output::Union{Num,Nothing}=nothing,
        solvetype::Symbol=:DISCRETE,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_name = Symbolics.tosymbol(param, escape=false)
        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end

        uhfunc = UHFunction(uhtype)

        #* Setup the name information of the hydroroutement
        meta = HydroMeta(inputs=[input_name], outputs=[output_name], params=[param_name], name=Symbol(output_name, :_uh_flux))

        return new{solvetype}(
            [input],
            [output],
            [param],
            uhfunc,
            meta
        )
    end
end


"""
    (flux::UnitHydroFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the unit hydrograph flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step (not supported, will throw an error).
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
- For matrix input: A matrix where each column is the result of applying the unit hydrograph to the corresponding input column.
- For array input: A array of routed outputs, with dimensions (output_var_names, node_names, ts_len).

# Notes
- The behavior differs based on the `solvetype` specified during the `UnitHydroFlux` construction:
  - `:unithydro1` uses a discrete problem solver approach.
  - `:unithydro2` uses a sparse matrix convolution approach.
- Vector input is not supported and will throw an error.
"""

(::UnitHydroFlux)(::Vector, ::ComponentVector; kwargs...) = @error "UnitHydroFlux is not support for single timepoint"

function (flux::UnitHydroFlux{:DISCRETE})(input::Matrix, pas::ComponentVector; kwargs...)
    solver = get(kwargs, :solver, DiscreteSolver())
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    input_vec = input[1, :]
    #* convert the lagflux to a discrete problem
    lag_prob!(du, u, p, t) = begin
        du[:] = input_vec[Int(t)] .* p[:weight] .+ [diff(u); -u[end]]
        nothing
    end
    #* prepare the initial states
    lag = pas[:params][get_param_names(flux)[1]]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))
    initstates = input_vec[1] .* uh_weight ./ sum(uh_weight)
    #* solve the problem
    sol = solver(lag_prob!, ComponentVector(weight=uh_weight ./ sum(uh_weight)), initstates, timeidx)
    reshape(sol[1, :], 1, length(input_vec))
end

function (flux::UnitHydroFlux{:SPARSE})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    lag = pas[:params][get_param_names(flux)[1]]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))
    #* the weight of the unit hydrograph is normalized by the sum of the weights
    uh_result = [-(i - 1) => uh_wi .* input_vec ./ sum(uh_weight) for (i, uh_wi) in enumerate(uh_weight)]
    #* construct the sparse matrix
    uh_sparse_matrix = spdiagm(uh_result...)
    #* sum the matrix
    sum_route = sum(uh_sparse_matrix, dims=2)[1:end-length(uh_weight)+1]
    reshape(sum_route, 1, length(input_vec))
end

# todo: 卷积计算的结果与前两个计算结果不太一致
function (flux::UnitHydroFlux{:INTEGRAL})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    itp_method = get(kwargs, :interp, LinearInterpolation)
    itp = itp_method(input_vec, collect(1:length(input_vec)), extrapolate=true)
    #* construct the unit hydrograph function based on the interpolation method and parameter
    lag = pas[:params][get_param_names(flux)[1]]
    tmax = get_uh_tmax(flux.uhfunc, lag)
    uh_sum = solve(IntegralProblem(flux.uhfunc, (0, tmax), lag), QuadGKJL()).u
    uh_itg_func = (x, p) -> flux.uhfunc(x, lag) * itp(p - x) / uh_sum
    #* solve the integral problem
    prob = IntegralProblem(uh_itg_func, (0, tmax), 1.0)
    routed_result = map(1:length(input_vec)) do t
        prob = remake(prob, p=t)
        sol = solve(prob, QuadGKJL())
        sol.u
    end
    reshape(routed_result, 1, length(input_vec))
end

# function (uh::AbstractUnitHydroFlux)(input::Array, pas::ComponentVector; kwargs...)
#     #* array dims: (variable dim, num of node, sequence length)
#     #* Extract the initial state of the parameters and routement in the pas variable
#     ptypes = get(kwargs, :ptypes, collect(keys(pas[:params])))
#     pytype_params = [pas[:params][ptype] for ptype in ptypes]

#     sols = map(eachindex(ptypes)) do (idx)
#         tmp_pas = ComponentVector(params=pytype_params[idx])
#         node_sols = reduce(hcat, uh(input[:, idx, :], tmp_pas))
#         node_sols
#     end
#     sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
#     return permutedims(sol_arr, (1, 3, 2))
# end