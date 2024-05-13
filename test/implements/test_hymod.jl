# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using ComponentArrays
using CairoMakie: Axis
using Interpolations

# test gr4j model
include("../../src/LumpedHydro.jl")

file_path = "data/hymod/LA011201_forcings.csv"
data = CSV.File(file_path);
df = DataFrame(data);
# time = 1:10000
prcp_vec = df[!, "precip"]
et_vec = df[!, "pet"]

# build model
s_max, b, a, kf, ks = 10, 5, 0.5, 0.5, 0.5
unit_params = ComponentVector(Smax=s_max, b=b, a=a, kf=kf, ks=ks)
init_states = ComponentVector(soilwater=0.0, fr1=1.0, fr2=1.0,
    fr3=1.0, sr=1.0
)
ps = ComponentVector(hymod=(params=unit_params, initstates=init_states, weight=1.0))

model = LumpedHydro.HyMOD.Node(name=:hymod, mtk=true)

input = (hymod=(prcp=prcp_vec, pet=et_vec,time=1:length(et_vec)),)
output = model(input, ps)
# output_df = DataFrames(output)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, 1:length(et_vec), output[:flow], color=:blue)
fig