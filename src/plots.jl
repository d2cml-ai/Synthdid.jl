
function synthdid_plot(estimates::synthdid_est1; treated_name="Treated", control_name="Synthetic control",
  treated_color="#043e7c", control_color="#d8450a",
  year_unit_trayectory=nothing,
  facet=nothing, lambda_comparable=!isnothing(facet), overlay=0,
  lambda_plot_scale=3, line_width=1, guide_linetype=:dash, point_size=3,
  diagram_alpha=0.95, trayectory_linewidth=2,
  se_method="jackknife")


  esti_mate = estimates.estimate

  expand_grid = collect(Base.product(1:length(esti_mate), 1:length(overlay)))
  s_estimate = [r[1] for r in expand_grid]
  s_overlay = [r[2] for r in expand_grid]


  vec_s = size(s_estimate, 2) * size(s_estimate, 1)

  s1_estimate = reshape(s_estimate, vec_s, 1)[:, 1]
  s1_overlay = reshape(s_overlay, vec_s, 1)[:, 1]
  grid = DataFrame(estimate=s1_estimate, overlay=s1_overlay)

  row = 1

  est = esti_mate[grid[row, "estimate"]]

  over = try
    overlay[grid[row, "overlay"]]
  catch
    overlay[row]
  end


  #TODO: implement vcov se
  se = if isnothing(se_method)
    return NaN
  else
    sqrt(vcov_synthdid_estimate(estimates, method=se_method))
  end


  setup, weights, N0, T0 = estimates.setup, estimates.weight, estimates.N0, estimates.T0
  Y = setup["Y"] .- contract3(setup["X"], weights["beta"])
  N1, T1 = size(Y) .- (N0, T0)

  lambda_synth = vcat(weights["lambda"], fill(0, T1))
  lambda_target = vcat(fill(0, T0), fill(1 / T1, T1))

  omega_synth = vcat(weights["omega"], fill(0, N1))
  omega_target = vcat(fill(0, N0), fill(1 / N1, N1))

  over = try
    over = estimates.overlat
  catch
    over = over
  end

  is_sc = all(weights["lambda"] .== 0) || over == 1

  intercept_offset = over * (omega_target .- omega_synth)' * Y * lambda_synth

  obs_trayectory = omega_target' * Y
  syn_trayectory = (omega_synth' * Y) .+ intercept_offset

  treated_post = omega_target' * Y * lambda_target
  treated_pre = omega_target' * Y * lambda_synth

  control_post = omega_synth' * Y * lambda_target + intercept_offset
  control_pre = omega_synth' * Y * lambda_synth + intercept_offset

  sdid_post = control_post + treated_pre - control_pre

  if isnothing(year_unit_trayectory)
    year_unit = collect(1:length(obs_trayectory))
  else
    year_unit = year_unit_trayectory
  end
  time = year_unit

  if length(time) == 0 || !all(isfinite.(time))
    time = 1:(T0+T1)
  end

  pre_time = lambda_synth' * time
  post_time = lambda_target' * time


  # treated_name = "Treated"
  # treated_color = "#d8450a"
  # control_name = "Synthetic control"
  # control_color = "#043e7c"

  lines = DataFrame(
    y=hcat(obs_trayectory, syn_trayectory)'[:, 1],
    x=repeat(time, 2),
    label=vcat(
      fill(treated_name, length(obs_trayectory)),
      fill(control_name, length(syn_trayectory))
    )
  )

  points = DataFrame(
    x=[post_time, post_time],
    y=[treated_post, sdid_post],
    label=[treated_name, control_name],
    color=[treated_color, control_color]
  )

  did_points = DataFrame(
    x=[pre_time, pre_time, post_time, post_time],
    y=[treated_pre, control_pre, control_post, treated_post],
    label=[treated_name, control_name, control_name, treated_name],
    color=[treated_color, control_color, control_color, treated_color]
  )

  hallucinated_segments = DataFrame(
    x=[pre_time, post_time],
    y=[treated_pre, sdid_post]
  )

  guide_segment = DataFrame(
    x=[pre_time, post_time, pre_time, post_time],
    y=[control_pre, control_post, treated_pre, sdid_post],
    label=["pre", "post", "pre", "post"]
  )

  arrow_df = DataFrame(
    x=[post_time, post_time],
    y=[sdid_post, treated_post],
    color=control_color
  )

  if lambda_comparable
    height = (maximum(obs_trayectory) - min(obs_trayectory)) / lambda_plot_scale
    bottom = minimum(obs_trayectory) .- height
    ribbons = DataFrame(x=time[1:T0], y=bottom .+ height .* lambda_synth[1:T0])
  else
    eval_algo = vcat(obs_trayectory', syn_trayectory')
    height = (
      maximum(eval_algo) -
      minimum(eval_algo)
    ) / lambda_plot_scale
    bottom = minimum(eval_algo) - height
    ribbons = DataFrame(
      x=time[1:T0],
      y=bottom .+ height .* lambda_synth[1:T0] ./ maximum(lambda_synth)
    )
  end

  T0s = try
    time[est.T0s]
  catch
    time[T0]
  end


  p = plot()
  plot!(time, obs_trayectory', label="Treated", color=treated_color, linewidth=trayectory_linewidth)
  plot!(time, syn_trayectory', label="Synthetic control", color=control_color, linewidth=trayectory_linewidth)
  scatter!(points.x, points.y, label="", color=points.color, ms=point_size, alpha=diagram_alpha)
  scatter!(did_points.x, did_points.y, color=did_points.color, label="", ms=point_size, alpha=diagram_alpha)
  plot!(hallucinated_segments.x, hallucinated_segments.y, color="black", label="", line=guide_linetype, linewidth=line_width)
  for i in [treated_name, control_name]
    df_algo = filter(row -> row.label == i, did_points)
    plot!(df_algo.x, df_algo.y, color=df_algo.color, label="")
  end
  for i in ["pre", "post"]
    df_algo = filter(row -> row.label == i, guide_segment)
    plot!(df_algo.x, df_algo.y, label="", color="black", line=guide_linetype, linewidth=line_width)
  end
  plot!()
  plot!(arrow_df.x, arrow_df.y, color=control_color, arrow=(0.5, 0.1), line_width=0.1, label="")
  plot!(ribbons.x, ribbons.y, fillrange=bottom, color=control_color, label="")
  plot!(ribbons.x, ribbons.y, color=treated_color, label="", linewidth=line_width)
  vline!([T0s], color="black", alpha=0.4, label="")
  # xlabel!("Time")
  # ylabel!("Trayectory")

  plot_description = Dict(
    "lines" => lines,
    "points" => points,
    "did_points" => did_points,
    "hallucinated_segments" => hallucinated_segments,
    "guide_segment" => guide_segment,
    "arrow_df" => arrow_df,
    "ribbons" => ribbons,
    "T0s" => T0s,
    "plot" => p)
  plot(p)
  return plot_description

end

function synthdid_units_plot(estimate::synthdid_est1; se_method::String="placebo", negligible_alpha::Float64=0.3, negligible_threshold::Float64=0.001, x_ticks=nothing)
  est = estimate

  setup, weights, N0, T0 = est.setup, est.weight, est.N0, est.T0

  Y = setup["Y"] - contract3(setup["X"], weights["beta"])

  N1, T1 = size(Y) .- (N0, T0)

  lambda_pre = vcat(weights["lambda"], fill(0, T1))
  lambda_post = vcat(zeros(T0), fill(1 / T1, T1))

  omega_control = vcat(weights["omega"], zeros(N1))
  omega_treat = vcat(zeros(N0), fill(1 / N1, N1))

  difs = omega_treat' * Y * (lambda_post .- lambda_pre) .- Y[1:N0, :] * (lambda_post .- lambda_pre)


  se = if isnothing(se_method)
    NaN
  else
    sqrt(vcov_synthdid_estimate(est, method=se_method))
  end

  include_units = if isnothing(x_ticks)
    1:N0
  else
    1:N0
  end

  est_imate = est.estimate

  plot_data = DataFrame(
    y=difs[include_units],
    x=include_units,
    weights=omega_control[include_units],
    estimate=est_imate,
    se=se
  )
  plot_data1 = filter(x -> x.weights > negligible_threshold, plot_data)
  plot_data2 = filter(x -> x.weights <= negligible_threshold, plot_data)

  p = plot()

  scatter!(plot_data1.x, plot_data1.y, ms=plot_data1.weights * 120, color="black", label="")
  scatter!(plot_data2.x, plot_data2.y, ms=plot_data2.weights * 120, color="black", label="", alpha=negligible_alpha,)
  hline!([est_imate], color="#000000", label="")
  if !isnan(se)
    hline!([est_imate - 1.96 * se], color="#000000", line=:dash, label="")
    hline!([est_imate + 1.96 * se], color="#000000", line=:dash, label="")
  end

  if !isnothing(x_ticks)
    xticks!(plot_data.x, x_ticks[include_units], rotation=90)
  end
  return p
end

# function synthdid_placebo_plot(estimate::synthdid_est1)
# end

function synthdid_rmse_plot(estimate::synthdid_est1)
  rmse = sqrt.(estimate.weight["vals"])

  data = DataFrame(rmse=rmse, x=1:length(rmse))
  p = plot()
  plot!(data.x, data.rmse, label="estimate")
  xlabel!("iteration")
  ylabel!("rmse")
  print(plot(p))
  return Dict("data" => data, "plot" => p)
end
