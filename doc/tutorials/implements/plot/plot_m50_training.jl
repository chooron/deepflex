using JLD2
using Plots
using Statistics
using Measures

m50_result = load("doc/tutorials/implements/save/m50_opt.jld2")
opt_loss_df = m50_result["loss_df"]
# Calculate moving average
function moving_average(data, window_size)
    output = similar(data)
    for i in 1:length(data)
        start_idx = max(1, i - window_size + 1)
        output[i] = mean(data[start_idx:i])
    end
    return output
end

# Create plot with customizations
p = plot(opt_loss_df.loss,
    label="Original Loss",
    color=:gray,
    alpha=0.4,
    xlabel="Iteration",
    ylabel="MSE",
    legendfont=font("times", 10),
    tickfont=font("times", 10),
    guidefont=font("times", 12),
    xticks=0:10:length(opt_loss_df.loss),
    dpi=300)

# Add moving average
window_size = 10  # Adjust this value to change the smoothing level
ma_loss = moving_average(opt_loss_df.loss, window_size)
plot!(ma_loss,
    label="Moving Average (n=$window_size)",
    color=:red,
    linewidth=2,
    alpha=0.8)
# savefig(p, "doc/tutorials/implements/plot/figures/m50_training_0.001.png")