
function plot_outcomes(res::Dict; save::Bool = false, dir::String = "", format::String = "png")

  if !(format in ["png", "pdf", "ps", "svg"])
    throw(ArgumentError("keyword argument `format` must be \"png\", \"pdf\", \"ps\", or \"svg\", got $format"))
  end
  
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
    # _intercept = (start_w - omega_hat)' * Y_pre_c * lambda_hat;
    Y_sdid_traj = omega_hat' * Y_c # .+ _intercept;

    # plot parameters for scaling/offsets
    plot_y_min = minimum([Y_sdid_traj; Y_t])
    plot_y_max = maximum([Y_sdid_traj; Y_t])
    plot_height = plot_y_max - plot_y_min

    p = Plots.plot(t_span, Y_sdid_traj', label = "Control", ls = :dash)
    plot!(t_span, Y_t', label = "Treatment")
    plot!(t_span[1:T0], lambda_hat .* plot_height' / 3 .+ plot_y_min, label = "", fillrange = plot_y_min, lw = 0)
    vline!([year], label = "")
    xlabel!("Year")
    title!("Adoption: " * year_str)
    plots[year_str] = p
    if save
      savefig(p, "$dir/outcomes_$year_str.$format")
    end
  end
  return plots
end

function plot_outcomes(
  data, Y_col::Union{String, Symbol}, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol}; 
  save::Bool = false, dir::String = "", format::String = "png", kwargs...
)
  if !(format in ["png", "pdf", "ps", "svg"])
    throw(ArgumentError("keyword argument `format` must be \"png\", \"pdf\", \"ps\", or \"svg\", got $format"))
  end

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...);
  return plot_outcomes(res, save = save, dir = dir, format = format)
end

function plot_weights(res::Dict; save::Bool = false, dir::String = "", format::String = "png")

  if !(format in ["png", "pdf", "ps", "svg"])
    throw(ArgumentError("keyword argument `format` must be \"png\", \"pdf\", \"ps\", or \"svg\", got $format"))
  end

  tyears = res["tyears"]
  year_params = res["year_params"]
  plots = Dict()

  for year in tyears
    year_str = string(year)
    # if isnothing(units) units = res["weights"]["omega"][:, S_col]
    units = res["weights"]["omega"][:, :country]
    omega_hat = res["weights"]["omega"][:, year_str]
    lambda_hat = res["weights"]["lambda"][year_str]
    tau_hat = year_params[year_params.treat_year .== year, :tau][1]
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
    plot_ms[plot_ms .== 24 * maximum(omega_hat)] .= 4
    
    p = Plots.plot(
        plot_units, plot_difs, seriestype = :scatter, label = "", mc = plot_cols, xrotation = 90, 
        xticks = (1:N0, units), ms = plot_ms, shape = plot_shape, msw = plot_msw, titlefontsize = 10, 
        tickfontsize = 6, labelfontsize = 8
    )
    hline!([tau_hat], label = "")
    xlabel!("Group")
    ylabel!("Difference")
    title!("Adoption: " * year_str)

    plots[year_str] = p
    if save
      savefig(p, "$dir/weights_$year_str.$format")
    end
  end

  return plots
end

function plot_weights(
  data, Y_col::Union{String, Symbol}, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol}; 
  save::Bool = false, dir::String = "", format::String = "png", kwargs...
)

  if !(format in ["png", "pdf", "ps", "svg"])
    throw(ArgumentError("keyword argument `format` must be \"png\", \"pdf\", \"ps\", or \"svg\", got $format"))
  end

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...)
  
  return plot_weights(res, save = save, dir = dir, format = format)
end