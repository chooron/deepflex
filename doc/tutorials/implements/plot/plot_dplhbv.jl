using JLD2

default(fontfamily="times")
# Create plot with customizations
p = plot(1:length(qobs_vec), qobs_vec,
    label="observed",
    color=:black,
    alpha=1.0,
    xlabel="Time Index",
    ylabel="Flow (mm)",
    legendfont=font("times", 10),
    dpi=300,
    xticks=0:1000:length(result.q))

plot!(1:length(result.q), result.q,
    label="simulated-before-optimization",
    alpha=0.6,
    color=:orange)

plot!(1:length(re_result.q), re_result.q,
    label="simulated-after-optimization",
    alpha=0.6,
    color=RGB(0.0, 0.4, 1.0))

mse_value1, fhv_value1, nse_value1 = mse(qobs_vec, result.q), fhv(qobs_vec, result.q, 0.1), nse(qobs_vec, result.q)
mse_value2, fhv_value2, nse_value2 = mse(qobs_vec, re_result.q), fhv(qobs_vec, re_result.q, 0.1), nse(qobs_vec, re_result.q)

# Add performance metrics as annotations
annotate!(length(result.q) / 3, maximum(qobs_vec) * 0.9,  # Move first annotation to the left third
    text("Before Optimization:\nMSE: $(round(mse_value1, digits=3))\nFHV: $(round(fhv_value1, digits=3))%\nNSE: $(round(nse_value1, digits=3))",
        :right, 10, "times"))

annotate!(2 * length(re_result.q) / 3, maximum(qobs_vec) * 0.9,  # Move second annotation to the right third
    text("After Optimization:\nMSE: $(round(mse_value2, digits=3))\nFHV: $(round(fhv_value2, digits=3))%\nNSE: $(round(nse_value2, digits=3))",
        :right, 10, "times"))

savefig("doc/tutorials/implements/plot/figures/dplHBV_predict.png")