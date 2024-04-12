#* copy from https://github.com/sclaw/muskingum-cunge/blob/main/src/muskingumcunge/reach.py
using DataInterpolations, NonlinearSolve, ModelingToolkit

struct BaseReach <: AbstractReach
    name::Symbol
    # reach attribute including width length slope mannings_n
    attr::NamedTuple
    # reach feature curve
    curve::NamedTuple
    upstream::Symbol
    downstream::Symbol
end

struct MuskingumReach <: AbstractReach
    name::Symbol
    attr::NamedTuple
    params::NamedTuple
    upstream::Symbol
    downstream::Symbol
end

# 基于一维水动力学的河道计算
struct HydrodynReach <: AbstractReach

end
"""
河宽,河长,坡降,糙率,河深
"""
function BaseReach(name::Symbol; attr::NamedTuple, upstream::Symbol, downstream::Symbol, resolution=50)
    stage = range(1e-4, get(attr, :max_stage, 10.0), resolution)
    width = repeat([get(attr, :width, 100.0)], resolution)
    area = stage .* width
    hydraulic_radius = @.(area / (width + (2 * stage)))
    discharge = @.((1 / get(attr, :mannings_n, 0.1)) * area * (hydraulic_radius^(2 / 3)) * get(attr, :slope, 100.0)^0.5)
    dq_da = (discharge[2:end] .- discharge[1:end-1]) / (area[2:end] .- area[1:end-1])
    celerity = push!(dq_da, dq_da[end])
    merge!(attr, (discharge=discharge, celerity=celerity))

    curve = (
        width=LinearInterpolation(log.(width), log.(discharge), extrapolate=true),
        celerity=LinearInterpolation(log.(celerity), log.(discharge), extrapolate=true)
    )


    BaseReach(
        name,
        attr,
        curve,
        upstream,
        downstream
    )
end

function TrapezoidalReach(name::Symbol; attr::NamedTuple, upstream::Symbol, downstream::Symbol, resolution=50)
    stage = range(0.0, get(attr, :max_stage, 10.0), resolution)
    #* calculate top width and use it as reach width
    width = @.(get(attr, :bottem_width, 100.0) + stage * get(attr, :side_slope, 100.0))
    area = @.((stage * get(attr, :side_slope, 100.0) + 2 * get(attr, :bottem_width, 100.0)) / 2 * stage)
    wetted_perimeter = @.(((stage * get(attr, :side_slope, 100.0))^2 + stage^2)^0.5 + get(attr, :bottem_width, 100.0))
    hydraulic_radius = @.(area / wetted_perimeter)
    discharge = @.((1 / get(attr, :mannings_n, 0.1)) * area * (hydraulic_radius^(2 / 3)) * get(attr, :slope, 100.0)^0.5)
    dq_da = (discharge[2:end] .- discharge[1:end-1]) ./ (area[2:end] .- area[1:end-1])
    celerity = push!(dq_da, dq_da[end])
    merge!(attr, (discharge=discharge, celerity=celerity))

    curve = (
        width=LinearInterpolation(log.(width), log.(discharge), extrapolate=true),
        celerity=LinearInterpolation(log.(celerity), log.(discharge), extrapolate=true)
    )

    BaseReach(
        name,
        attr,
        curve,
        upstream,
        downstream
    )
end

function RectangleReach(name; attr, upstream::Symbol, downstream::Symbol, resolution=50)
    stage = range(1e-4, get(attr, :max_stage, 10.0), resolution)
    #* calculate top width and use it as reach width
    area = stage .* get(attr, :width, 100.0)
    alpha = @.((1 / get(attr, :mannings_n, 0.01)) * (get(attr, :slope, 100.0)^0.5) * (get(attr, :width, 100.0)^(-2 / 3)))
    discharge = alpha .* area .* (5 / 3)

    #* calculate muskingum params k and x
    c = alpha * (5 / 3) * (area^(5 / 3 - 1))
    c[1] = c[2]
    k = get(attr, :length, 1000.0) / c / 3600 # unit is hours 
    x = (1 / 2) - discharge / (2 * c * get(attr, :width, 100.0) * get(attr, :slope, 100.0) * get(attr, :length, 1000.0))

    MuskingumReach(
        name,
        attr,
        (k=k, x=x),
        upstream,
        downstream
    )
end

function calcu_muskingum_params(reach::BaseReach, reach_q::Number, dt::Number)
    #* 根据河道在t~t+1时段内的平均流量结合特征曲线计算调整后的参数
    log_reach_q = log(abs(reach_q))

    tmp_b = exp(reach.curve[:width](log_reach_q))
    tmp_c = exp(reach.curve[:celerity](log_reach_q))

    courant = tmp_c * dt * 60 * 60 / reach.attr.length
    reynold = reach_q / (reach.attr.slope * tmp_c * reach.attr.length * tmp_b)

    c0 = (-1 + courant + reynold) / (1 + courant + reynold)
    c1 = (1 + courant - reynold) / (1 + courant + reynold)
    c2 = (1 - courant + reynold) / (1 + courant + reynold)

    #* 获取计算后的参数
    c0, c1, c2
end

function (reach::BaseReach)(input::AbstractVector, dt::Number=1.0)
    output = zeros(eltype(input), size(input))
    output[1] = input[1]
    min_input = minimum(input)

    @assert maximum(input) < maximum(reach.attr.discharge)

    #* build nolinear problem
    function nolinear_func(u, p)
        inflow_t1, inflow_t2, outflow_t1 = p[1], p[2], p[3]
        c0, c1, c2 = calcu_muskingum_params(reach, (inflow_t1 + inflow_t2 + outflow_t1 + u[1]) / 4, dt)
        q_guess = c0 * inflow_t2 + c1 * inflow_t1 + c2 * outflow_t1
        u .- max(min_input, q_guess)
    end

    #* iter solve
    for i in 1:(length(input)-1)
        guess = [(input[i] + input[i+1] + output[i]) / 3]
        ps = [
            input[i],
            input[i+1],
            output[i],
        ]
        prob = NonlinearProblem(nolinear_func, guess, ps)

        sol = solve(prob, NewtonRaphson())
        output[i+1] = first(sol.u)
    end
    output
end

function (reach::MuskingumReach)(input::AbstractVector, dt::Number=1.0)
    output = zeros(eltype(input), size(input))
    output[1] = input[1]
    k, x = reach.params.k, reach.params.x

    c0 = ((dt / k) - (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c1 = ((dt / k) + (2 * x)) / ((2 * (1 - x)) + (dt / k))
    c2 = ((2 * (1 - x)) - (dt / k)) / ((2 * (1 - x)) + (dt / k))

    for i in 1:(length(input)-1)
        q_out = (c0 * input[i+1]) + (c1 * input[i]) + (c2 * outflows[i])
        output[i+1] = max(minimum(input), q_out)
    end
    output
end
