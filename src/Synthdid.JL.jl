module Synthdid

export data, estimate_dgp, simulate_dgp, synthdid_plot, synthdid_units_plot, synthdid_rmse_plot, fw_step, sc_weight_fw, summary_synth, synthdid_estimate, sc_estimate, did_estimate, synthdid_effect_curve, panel_matrices, random_low_rank, vcov_synthdid_estimate

using DataFrames, Plots, CSV, Statistics, Shuffle, Distributions

# Write your package code here.
include("data.jl")
include("utils.jl")
include("placebo_simulations.jl")
include("solver.jl")
include("synthdid.jl")
include("vcov.jl")
include("plots.jl")
include("summary.jl")
end
