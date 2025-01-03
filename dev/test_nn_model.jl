# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using Plots
using Lux
using StableRNGs

include("../src/HydroModels.jl")
HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
NeuralFlux = HydroModels.NeuralFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/m50.jl")