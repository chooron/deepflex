#* copy from https://github.com/sclaw/muskingum-cunge/blob/main/src/muskingumcunge/reach.py
using DataInterpolations
abstract type AbstractReach end
struct BaseReach <: AbstractReach
    width
    length
    slope
    mannings_n
    discharge
    celerity_curve
    width_curve
end

struct MuskingumReach <: AbstractReach
    k
    x
end

"""
河宽,河长,坡降,糙率,河深
"""
function BaseReach(width, length, slope, mannings_n, max_stage=10, stage_res=50)
    stage = range(1e-4, max_stage, stage_res)
    width = repeat([width], stage_res)
    area = stage .* width
    hydraulic_radius = @.(area / (width + (2 * stage)))
    discharge = @.((1 / mannings_n) * area * (hydraulic_radius^(2 / 3)) * slope^0.5)
    dq_da = (discharge[2:end] .- discharge[1:end-1]) / (area[2:end] .- area[1:end-1])
    celerity = push!(dq_da, dq_da[end])

    width_curve = LinearInterpolation(log.(width), log.(discharge), extrapolate=true)
    celerity_curve = LinearInterpolation(log.(celerity), log.(discharge), extrapolate=true)

    BaseReach(
        width,
        length,
        slope,
        mannings_n,
        discharge,
        width_curve,
        celerity_curve
    )
end

function TrapezoidalReach(bottem_width, length, slope, mannings_n, side_slope, max_stage=10, stage_res=50)
    stage = range(1e-4, max_stage, stage_res)
    #* calculate top width and use it as reach width
    width = @.(bottem_width + stage * side_slope)
    area = @.((stage * side_slope + 2 * bottem_width) / 2 * stage)
    wetted_perimeter = @.(((stage * side_slope)^2 + stage^2)^0.5 + bottem_width)
    hydraulic_radius = @.(area / wetted_perimeter)

    discharge = @.((1 / mannings_n) * area * (hydraulic_radius^(2 / 3)) * slope^0.5)
    dq_da = (discharge[2:end] .- discharge[1:end-1]) ./ (area[2:end] .- area[1:end-1])
    celerity = push!(dq_da, dq_da[end])

    width_curve = LinearInterpolation(log.(width), log.(discharge), extrapolate=true)
    celerity_curve = LinearInterpolation(log.(celerity), log.(discharge), extrapolate=true)

    BaseReach(
        width,
        length,
        slope,
        mannings_n,
        discharge,
        width_curve,
        celerity_curve
    )
end

function RectangleReach(width, length, slope, mannings_n, max_stage=10, stage_res=50)
    stage = range(1e-4, max_stage, stage_res)
    #* calculate top width and use it as reach width
    area = stage .* width
    alpha = @.((1 / mannings_n) * (slope^0.5) * (width^(-2 / 3)))
    discharge = alpha .* area .* (5 / 3)

    #* calculate muskingum params k and x
    c = alpha * (5 / 3) * (area^(5 / 3 - 1))
    c[1] = c[2]
    k = length / c / 3600 # unit is hours 
    x = (1 / 2) - discharge / (2 * c * width * slope * length)

    MuskingumReach(
        k,
        x
    )
end

function calcu_muskingum_params(reach::BaseReach, reach_q::Number, dt::Number)
    #* 根据河道在t~t+1时段内的平均流量结合特征曲线计算调整后的参数
    courant = exp(reach.celerity_curve(reach_q)) * dt * 60 * 60 / reach.length
    reynold = reach_q / (reach.slope * reach.celerity_curve(reach_q) * reach.length * reach.width_curve(reach_q))

    c0 = (-1 + courant + reynold) / (1 + courant + reynold)
    c1 = (1 + courant - reynold) / (1 + courant + reynold)
    c2 = (1 - courant + reynold) / (1 + courant + reynold)

    #* 获取计算后的参数
    c0, c1, c2
end

function (reach::BaseReach)(input::AbstractVector, dt::Number=1.0)
    output = zeros(eltype(input), size(input))
    output[1] = input[1]

    @assert maximum(input) < maximum(reach.discharge)

    #* build nolinear problem
    @variables outflow_t2
    @parameters inflow_t1, inflow_t2, outflow_t1

    eqs = [outflow_t2 ~ begin
        c0, c1, c2 = calcu_muskingum_params(reach, (inflow_t1 + inflow_t2 + outflow_t1 + outflow_t2) / 4, dt)
        q_guess = c0 * inflow_t2 + c1 * inflow_t1 + c2 * outflow_t1
        max(minimum(input), q_guess)
    end]

    @mtkbuild ns = NonlinearSystem(eqs, [outflow_t2], [inflow_t1, inflow_t2, outflow_t1])

    #* iter solve
    for i in 1:(length(input)-1)
        guess = [outflow_t2 => (input[i] + input[i+1] + output[i]) / 3]
        ps = [
            inflow_t1 => input[i],
            inflow_t2 => input[i+1],
            outflow_t1 => output[i],
        ]
        prob = NonlinearProblem(ns, guess, ps)
        sol = solve(prob, NewtonRaphson())
        println(sol)
        output[i+1] = first(sol.u)
    end
    output
end

function (reach::MuskingumReach)(input::AbstractVector, dt::Number=1.0)
    output = zeros(eltype(input), size(input))
    output[1] = input[1]

    c0 = ((dt / reach.k) - (2 * reach.x)) / ((2 * (1 - reach.x)) + (dt / reach.k))
    c1 = ((dt / reach.k) + (2 * reach.x)) / ((2 * (1 - reach.x)) + (dt / reach.k))
    c2 = ((2 * (1 - reach.x)) - (dt / reach.k)) / ((2 * (1 - reach.x)) + (dt / reach.k))

    for i in 1:(length(input)-1)
        q_out = (c0 * input[i+1]) + (c1 * input[i]) + (c2 * outflows[i])
        q_out = max(minimum(input), q_out)
        output[i+1] = q_out
    end
end
