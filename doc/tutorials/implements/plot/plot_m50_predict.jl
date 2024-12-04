using JLD2
using Plots
using Measures
using Statistics
using Printf
# include("tutorials/implements/loss_functions.jl")

result = load("tutorials/implements/save/m50_result.jld2")
obs = result["obs"]
exphydro = result["exphydro"]
qnn = result["qnn"]
m50 = result["m50"]

mse_value1, fhv_value1, nse_value1 = mse(obs, exphydro), fhv(obs, exphydro, 0.1), nse(obs, exphydro)
mse_value2, fhv_value2, nse_value2 = mse(obs, qnn), fhv(obs, qnn, 0.1), nse(obs, qnn)
mse_value3, fhv_value3, nse_value3 = mse(obs, m50), fhv(obs, m50, 0.1), nse(obs, m50)
criterion = [
    [mse_value1, fhv_value1, nse_value1],
    [mse_value2, fhv_value2, nse_value2],
    [mse_value3, fhv_value3, nse_value3]
]

# 创建三行一列的子图
function plot_predictions(result)
    # 设置全局字体
    default(fontfamily="times", legendfontsize=14,
        tickfontsize=12, guidefontsize=14)

    p1 = plot(
        layout=(3, 1), size=(800, 1000), margin=0.2mm,
        bottom_margin=2mm, top_margin=2mm, left_margin=5mm, right_margin=5mm,
        link=:x, dpi=600
    )
    time = collect(1:10000)
    # 定义颜色和透明度
    colors = [RGB(1.0, 0.4, 0.0), :green, RGB(0.0, 0.4, 1.0)]  # exphydro, qnn, m50的颜色

    # 三个子图的数据和透明度设置
    alphas = [
        [0.8, 0.3, 0.3],  # 第一个子图：exphydro不透明
        [0.3, 0.8, 0.3],  # 第二个子图：qnn不透明
        [0.3, 0.3, 0.8]   # 第三个子图：m50不透明
    ]

    labels = ["(a)", "(b)", "(c)"]
    for i in 1:3
        # 首先绘制观测值（在最底层）
        plot!(p1[i], time, result["obs"],
            label="Observed", color=:black, alpha=1.0,
            subplot=i)

        # 按照透明度顺序绘制模型预测值
        models = ["exphydro", "qnn", "m50"]
        model_names = ["ExpHydro", "NNq", "M50"]  # 更好的显示名称
        for (j, (model, name)) in enumerate(zip(models, model_names))
            plot!(p1[i], time, result[model],
                label=name, color=colors[j],
                alpha=alphas[i][j], linewidth=2,
                subplot=i)
        end

        # 设置子图标签和标注
        plot!(p1[i], xlabel=(i == 3 ? "Time Index" : ""),
            ylabel="Flow(mm)",
            title="",
            annotations=(50, maximum(result["obs"]) * 1.1,
                text(labels[i], :left, 14, "times")))

        # 添加评价指标文本
        criterion_text = "MSE: $(round(criterion[i][1], digits=3))\nFHV: $(round(criterion[i][2], digits=3))%\nNSE: $(round(criterion[i][3], digits=3))"
        annotate!(p1[i], 5000, maximum(result["obs"]) * 0.95,
            text(criterion_text, :left, 14, "times"))
    end

    return p1
end

# 绘制图形
p = plot_predictions(result)
savefig(p, "tutorials/implements/save/m50_predict.png")