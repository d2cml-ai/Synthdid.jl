
function plot_outcomes(
  data, Y_col::Union{String, Symbol}, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol}; kwargs...
)

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...);
  t_span = res["t_span"]; 
  tyears = res["tyears"];
  N0 = size(res["weights"]["omega"], 1)
  plots = Dict()

  for year in tyears
    year_str = string(year);
    T0 = res["year_params"][res["year_params"].treat_year .== year, "T0"][1];
    T0 = Int(T0);

    # get weights
    omega_hat = res["weights"]["omega"][:, year_str];
    lambda_hat = res["weights"]["lambda"][year_str];

    # create matrices to calculate trajectory
    Y_year = copy(res["Y"][year_str]);
    Y_pre_c = Y_year[1:N0, 1:T0];
    Y_post_c = Y_year[1:N0, T0 + 1:end];
    Y_pre_t = mean(Y_year[N0 + 1:end, 1:T0], dims = 1);
    Y_post_t = mean(Y_year[N0 + 1:end, T0 + 1:end], dims = 1);
    Y_t = [Y_pre_t Y_post_t];
    Y_c = [Y_pre_c Y_post_c]

    # calculate synth outcome trajectory
    n_features = size(Y_c, 1);
    start_w = fill(1/n_features, n_features);
    omega_hat = res["weights"]["omega"][:, year_str];
    lambda_hat = res["weights"]["lambda"][year_str];
    _intercept = (start_w - omega_hat)' * Y_pre_c * lambda_hat;
    Y_sdid_traj = omega_hat' * Y_c # .+ _intercept;

    # plot parameters for scaling/offsets
    plot_y_min = minimum([Y_sdid_traj; Y_t])
    plot_y_max = maximum([Y_sdid_traj; Y_t])
    plot_height = plot_y_max - plot_y_min

    p = plot(t_span, Y_sdid_traj', label = "Control", ls = :dash)
    plot!(t_span, Y_t', label = "Treatment")
    plot!(t_span[1:T0], lambda_hat .* plot_height' / 3 .+ plot_y_min, label = "", fillrange = plot_y_min, lw = 0)
    vline!([year], label = "")
    xlabel!("Year")
    title!("Adoption: " * year_str)
    plots[year_str] = p
  end
  return plots
end

function plot_weights(
  data, Y_col::Union{String, Symbol}, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol}; kwargs...
)

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...)
  tyears = res["tyears"]
  plots = Dict()

  for year in tyears
    year_str = string(year)
    # if isnothing(units) units = res["weights"]["omega"][:, S_col]
    units = res["weights"]["omega"][:, :country]
    omega_hat = res["weights"]["omega"][:, year_str]
    lambda_hat = res["weights"]["lambda"][year_str]
    N0 = size(omega_hat, 1)
    T0 = size(lambda_hat, 1)
    Y = copy(res["Y"][year_str])
    N, T = size(Y)
    N1, T1 = (N, T) .- (N0, T0)

    lambda_pre = [lambda_hat; zeros(T1)]
    lambda_post = [zeros(T0); fill(1/T1, T1)]
    omega_control = [omega_hat; zeros(N1)]
    omega_treat = [zeros(N0); fill(1/N1, N1)]
    difs = omega_treat' * Y * (lambda_post - lambda_pre) .- Y[1:N0, :] * (lambda_post - lambda_pre)

    plot_units = reshape(units, (1, N0))
    plot_difs = reshape(difs, (1, N0))
    plot_ms = reshape(omega_hat, (1, N0)) / maximum(omega_hat)
    plot_shape = @. ifelse(plot_ms == 0, :xcross, :circle)
    plot_cols = @. ifelse(plot_ms == 0, :red, :dodgerblue4)
    plot_msw = @. ifelse(plot_ms == 0, 1, 0)
    plot_ms = (plot_ms .+ 4 * maximum(omega_hat)) * 6
    plot_ms[plot_ms .== 40 * maximum(omega_hat)] .= 4
    
    p = plot(
        plot_units, plot_difs, seriestype = :scatter, label = "", mc = plot_cols, xrotation = 90, 
        xticks = (1:N0, units), ms = plot_ms, shape = plot_shape, msw = plot_msw, titlefontsize = 10, 
        tickfontsize = 6, labelfontsize = 8
    )
    hline!([res["att"]], label = "")
    xlabel!("Group")
    ylabel!("Difference")
    title!("Adoption: " * year_str)

    plots[year_str] = p
  end

  return plots
end

# function synthdid_units_plot(estimate::synthdid_est1; se_method::String="placebo", negligible_alpha::Float64=0.3, negligible_threshold::Float64=0.001, x_ticks=nothing)
#   est = estimate

#   setup, weights, N0, T0 = est.setup, est.weight, est.N0, est.T0

#   Y = setup["Y"] - contract3(setup["X"], weights["beta"])

#   N1, T1 = size(Y) .- (N0, T0)

#   lambda_pre = vcat(weights["lambda"], fill(0, T1))
#   lambda_post = vcat(zeros(T0), fill(1 / T1, T1))

#   omega_control = vcat(weights["omega"], zeros(N1))
#   omega_treat = vcat(zeros(N0), fill(1 / N1, N1))

#   difs = omega_treat' * Y * (lambda_post .- lambda_pre) .- Y[1:N0, :] * (lambda_post .- lambda_pre)


#   se = if isnothing(se_method)
#     NaN
#   else
#     sqrt(vcov_synthdid_estimate(est, method=se_method))
#   end

#   include_units = if isnothing(x_ticks)
#     1:N0
#   else
#     1:N0
#   end

#   est_imate = est.estimate

#   plot_data = DataFrame(
#     y=difs[include_units],
#     x=include_units,
#     weights=omega_control[include_units],
#     estimate=est_imate,
#     se=se
#   )
#   plot_data1 = filter(x -> x.weights > negligible_threshold, plot_data)
#   plot_data2 = filter(x -> x.weights <= negligible_threshold, plot_data)

#   p = plot()

#   scatter!(plot_data1.x, plot_data1.y, ms=plot_data1.weights * 120, color="black", label="")
#   scatter!(plot_data2.x, plot_data2.y, ms=plot_data2.weights * 120, color="black", label="", alpha=negligible_alpha,)
#   hline!([est_imate], color="#000000", label="")
#   if !isnan(se)
#     hline!([est_imate - 1.96 * se], color="#000000", line=:dash, label="")
#     hline!([est_imate + 1.96 * se], color="#000000", line=:dash, label="")
#   end

#   if !isnothing(x_ticks)
#     xticks!(plot_data.x, x_ticks[include_units], rotation=90)
#   end
#   return p
# end

# # function synthdid_placebo_plot(estimate::synthdid_est1)
# # end

# function synthdid_rmse_plot(estimate::synthdid_est1)
#   rmse = sqrt.(estimate.weight["vals"])

#   data = DataFrame(rmse=rmse, x=1:length(rmse))
#   p = plot()
#   plot!(data.x, data.rmse, label="estimate")
#   xlabel!("iteration")
#   ylabel!("rmse")
#   print(plot(p))
#   return Dict("data" => data, "plot" => p)
# end
