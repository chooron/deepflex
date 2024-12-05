"""
 mean squared error
"""
mse(obs,sim) = mean((obs .- sim).^2)

""" 
 nash sutcliffe efficient 
"""
nse(obs,sim) = 1 - sum((obs .- sim).^2)/sum((obs .- mean(obs)).^2)

"""
    fhv(obs, sim, h=0.02)

Calculate the high flow volume error (FHV) metric.
h: fraction of highest flows to consider (default is 0.02 or top 2%)
"""
function fhv(obs, sim, h=0.02)
    # Sort discharges in descending order
    obs_sorted = sort(obs, rev=true)
    sim_sorted = sort(sim, rev=true)
    
    # Subset data to only top h flow values
    n = round(Int, h * length(obs))
    obs_high = obs_sorted[1:n]
    sim_high = sim_sorted[1:n]
    
    # Calculate FHV
    return sum(sim_high - obs_high) / sum(obs_high) *100
end