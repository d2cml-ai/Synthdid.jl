
mutable struct SynthDID

  data::DataFrame
  Y_col::Symbol
  S_col::Symbol
  T_col::Symbol
  D_col::Symbol
  proc_data::DataFrame
  covariates::Union{Nothing, Vector}
  tyears::Vector
  t_span::Vector
  omega_hat::DataFrame
  lambda_hat::Dict
  beta_hat::Union{Nothing, DataFrame}
  ATT::Float64
  se::Float64
  t_stat::Float64
  p_val::Float64
  Y::Dict
  X::Union{Dict, Nothing}
  N::Int
  T::Int
  N0::Int
  T0::DataFrame
end

function SynthDID(data, Y_col, S_col, T_col, D_col; se_method::String = "placebo", se_reps::Int = 50, kwargs...)

  # ensure se_method is a valid option
  if !(se_method in ["jackknife", "bootstrap", "placebo"])
    throw(ArgumentError("se_method must be 'jackknife', 'bootstrap', or 'placebo', got $se_method"))
  end

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...)
  proc_data = res["proc_data"]
  covariates = res["covariates"]
  tyears = res["tyears"]
  t_span = res["t_span"]
  omega_hat = res["weights"]["omega"]
  lambda_hat = res["weights"]["lambda"]
  beta_hat = res["weights"]["beta"]
  ATT = res["att"]
  Y = res["Y"]
  X = res["X"]
  N, T = res["N"], res["T"]
  N0, T0 = res["N0"], res["T0"]

  if se_method == "jackknife"
    se = jackknife_se(data, Y_col, S_col, T_col, D_col, att = ATT, omega_hat = omega_hat, lambda_hat = lambda_hat, kwargs...)
  elseif se_method == "bootstrap"
    se = bootstrap_se(data, Y_col, S_col, T_col, D_col, n_boot = se_reps, kwargs...)
  else
    se = placebo_se(data, Y_col, S_col, T_col, D_col, n_reps = se_reps, kwargs...)
  end

  t_stat = ATT / se
  p_val = 2 * (1 - cdf(Normal(0, 1), abs(t_stat)))

  obj = SynthDID(data, Y_col, S_col, T_col, D_col, proc_data, covariates, tyears, t_span, omega_hat, lambda_hat, beta_hat, ATT, se, t_stat, p_val, Y, X, N, T, N0, T0)
  return obj
end

  

  
