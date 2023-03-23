module Synthdid

export fw_step, sc_weight_fw, sc_weight_covariates, sdid, california_prop99, quota

using DataFrames, Plots, CSV, Statistics, Shuffle, Distributions

# Write your package code here.
include("data.jl")
include("utils.jl")
# include("placebo_simulations.jl")
include("solver.jl")
include("main.jl")
# include("vcov.jl")
# include("plots.jl")
# include("summary.jl")
end