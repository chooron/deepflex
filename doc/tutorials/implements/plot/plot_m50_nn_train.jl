using JLD2
using Plots
using Statistics
using Measures

# Load the optimization results
m50_nn_opt_result = load("doc/tutorials/implements/save/m50_nn_opt.jld2")

epnn_loss_df = m50_nn_opt_result["epnn_loss_df"]
qnn_loss_df = m50_nn_opt_result["qnn_loss_df"]

default(fontfamily="times")
# Calculate moving average
function moving_average(data, window_size)
    output = similar(data)
    for i in 1:length(data)
        start_idx = max(1, i - window_size + 1)
        output[i] = mean(data[start_idx:i])
    end
    return output
end


function plot_losses(epnn_loss_df, qnn_loss_df)
    # Create a figure with two subplots
    p = plot(layout=(1, 2), size=(900, 300), dpi=300,
        left_margin=5mm,    # 增加左边距
        bottom_margin=5mm,  # 增加底部边距
        right_margin=5mm)

    # Plot EPNN loss (subplot a)
    window_size = 10
    plot!(p[1], epnn_loss_df.loss,
        label="Original Loss",
        color=:gray,
        alpha=0.4,
        xlabel="Iteration",
        ylabel="MSE",
        title="(a) NNep Training Loss",
        legendfont=font("times", 10),
        tickfont=font("times", 10),
        guidefont=font("times", 12))

    ma_loss = moving_average(epnn_loss_df.loss, window_size)
    plot!(p[1], ma_loss,
        label="Moving Average (n=$window_size)",
        color=:red,
        linewidth=2,
        alpha=0.8)

    # Plot QNN loss (subplot b)
    plot!(p[2], qnn_loss_df.loss,
        label="Original Loss",
        color=:gray,
        alpha=0.4,
        xlabel="Iteration",
        ylabel="MSE",
        title="(b) NNq Training Loss",
        legendfont=font("times", 10),
        tickfont=font("times", 10),
        guidefont=font("times", 12))

    ma_loss = moving_average(qnn_loss_df.loss, window_size)
    plot!(p[2], ma_loss,
        label="Moving Average (n=$window_size)",
        color=:red,
        linewidth=2,
        alpha=0.8)

    # Save the figure
    savefig(p, "doc/tutorials/implements/plot/figures/m50_nn_training.png")
end

# Replace the original plot_loss calls with:
plot_losses(epnn_loss_df, qnn_loss_df)

