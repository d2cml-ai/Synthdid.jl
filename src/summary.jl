mutable struct summary_synthdid
  estimate
  se
  control
  periods
end


mutable struct summary_element
  N1
  T1
  N0
  T0
  effective
  data
end

function synthdid_controls(estimates; desc=true, weight_type="omega", panel)
  weight = estimates.weight[weight_type]
  eff = 1 / sum(weight .^ 2)
  n0 = estimates.N0
  t0 = estimates.T0
  Y = estimates.setup["Y"]
  n1, t1 = size(Y) .- (n0, t0)

  time = panel.time

  info_weight = DataFrame(weight=weight)
  if weight_type != "omega"
    info_weight.time = time[1:t0]
    sort!(info_weight, :weight, rev=desc)
    filter!(x -> x.weight != 0, info_weight)
  else
    info_weight = sort!(info_weight, :weight, rev=desc)[1:Int(floor(size(info_weight, 1) / 2)), :]
    filter!(x -> x.weight != 0, info_weight)
  end
  return summary_element(n1, t1, n0, t0, eff, info_weight)
end


function summary_synth(estimates::synthdid_est1, print_result=true; panel, print_all=false, round_digits=3)
  estimate = estimates.estimate
  # panel_matrices = panel_matrices
  se = vcov_synthdid_estimate(estimates)
  control = synthdid_controls(estimates, panel=panel)
  periods = synthdid_controls(estimates, weight_type="lambda", panel=panel)

  if print_result
    est = round(estimate, digits=round_digits)
    eff_n0 = round(control.effective, digits=round_digits)
    n0_eff = round(eff_n0 / control.N0, digits=round_digits)
    eff_t0 = round(periods.effective, digits=round_digits)
    t0_eff = round(eff_t0 / periods.T0, digits=round_digits)

    summary_result = "synthdid: $(est) +- $(1.96 * se). Effective N0/N0 = $(eff_n0)/$(control.N0)~$(n0_eff). Effective T0/T0 = $(eff_t0)/$(periods.T0) ~ $(t0_eff). N1,T1 = $(control.N1), $(periods.T1)."
    print(summary_result)
  end

  if print_all
    print(control.data)
    print(periods.data)
  end
  return summary_synthdid(estimate, se, control, periods)
end
