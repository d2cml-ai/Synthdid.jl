module Synthdid

export fw_step, sc_weight_fw, sc_weight_covariates, sdid, california_prop99, quota, jackknife_se, bootstrap_se, placebo_se, plot_outcomes, plot_weights, SynthDID

using DataFrames, Plots, CSV, Statistics, Distributions, Random

# Write your package code here.
include("main.jl")
include("data.jl")
include("utils.jl")
# include("placebo_simulations.jl")
include("solver.jl")
include("estimations.jl")
include("vcov.jl")
include("plots.jl")
# include("summary.jl")
end