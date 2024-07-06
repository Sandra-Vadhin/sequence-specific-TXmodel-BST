# packages -
using CSV
using DataFrames
using LinearAlgebra
using Plots
using Distributions
using DifferentialEquations
using GlobalSensitivity
using Tables
using NumericalIntegration
using Optim
using Colors
using JLD2
 using BSTModelKit
 using ForwardDiff, DiffResults
using Sundials
using DelimitedFiles
using SmoothingSplines
using Interpolations
# using Optimization
# using OptimizationOptimJL
# using OptimizationBBO


# paths -
_PATH_TO_MODEL = joinpath(pwd(),"model")
_PATH_TO_DATA = joinpath(pwd(),"data")
_PATH_TO_TMP = joinpath(pwd(),"tmp")
_PATH_TO_SRC = joinpath(pwd(),"src")
_PATH_TO_FIGS = joinpath(pwd(),"figs")
_PATH_TO_ACTUAL_ENSEMBLE = joinpath(pwd(),"actual_ensemble_s_system")

# code -

#include(joinpath(_PATH_TO_SRC, "Include.jl"))
include(joinpath(_PATH_TO_SRC,"Evaluate.jl"))
 include(joinpath(_PATH_TO_SRC,"Balance_Opt.jl"))
include(joinpath(_PATH_TO_SRC,"Kinetic.jl"))
include(joinpath(_PATH_TO_SRC,"Learn.jl"))
include(joinpath(_PATH_TO_SRC,"Compute.jl"))
# include("thrombin.jl")
include(joinpath(_PATH_TO_SRC,"Loss.jl"))