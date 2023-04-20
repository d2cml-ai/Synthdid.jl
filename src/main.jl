
mutable struct SynthDID

  data::DataFrame
  Y_col::Symbol
  S_col::Symbol
  T_col::Symbol
  D_col::Symbol
  proc_data::DataFrame
  covariates::Union{Nothing, Vector}
  cov_method::Union{String, Nothing}
  se_method::String
  tyears::Vector
  t_span::Vector
  omega_hat::DataFrame
  lambda_hat::Dict
  beta_hat::Union{Nothing, DataFrame}
  year_params::DataFrame
  ATT::Float64
  se::Float64
  t_stat::Float64
  p_val::Float64
  Y::Dict
  X::Union{Dict, Nothing, Matrix}
  N::Int
  T::Int
  N0::Int
  T0::DataFrame
  coef_table::DataFrame
end

function SynthDID(data, Y_col, S_col, T_col, D_col; se_method::String = "placebo", se_reps::Int = 50, kwargs...)

  # ensure se_method is a valid option
  if !(se_method in ["jackknife", "bootstrap", "placebo"])
    throw(ArgumentError("se_method must be 'jackknife', 'bootstrap', or 'placebo', got $se_method"))
  end

  Y_col = Symbol(Y_col)
  S_col = Symbol(S_col)
  T_col = Symbol(T_col)
  D_col = Symbol(D_col)

  res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...)
  proc_data = res["proc_data"]
  covariates = res["covariates"]
  cov_method = res["cov_method"]
  tyears = res["tyears"]
  t_span = res["t_span"]
  omega_hat = res["weights"]["omega"]
  lambda_hat = res["weights"]["lambda"]
  beta_hat = res["weights"]["beta"]
  year_params = res["year_params"]
  ATT = res["att"]
  Y = res["Y"]
  X = res["X"]
  N, T = res["N"], res["T"]
  N0, T0 = res["N0"], res["T0"]

  if se_method == "jackknife"
    se = jackknife_se(data, Y_col, S_col, T_col, D_col; att = ATT, omega_hat = omega_hat, lambda_hat = lambda_hat, kwargs...)
  elseif se_method == "bootstrap"
    se = bootstrap_se(data, Y_col, S_col, T_col, D_col; n_boot = se_reps, kwargs...)
  else
    se = placebo_se(data, Y_col, S_col, T_col, D_col; n_reps = se_reps, kwargs...)
  end

  t_stat = ATT / se
  p_val = 2 * (1 - cdf(Normal(0, 1), abs(t_stat)))

  table_cols = ["ATT", "Std. Err.", "t", "P>|t|", "[95% Conf.", "Interval]"]
  ci_inf = ATT - se * 1.96
  ci_sup = ATT + se * 1.96
  coef_table = DataFrame([[ATT], [se], [t_stat], [p_val], [ci_inf], [ci_sup]], table_cols)
  
  obj = SynthDID(
    data, Y_col, S_col, T_col, D_col, proc_data, covariates, cov_method, se_method, tyears, t_span, omega_hat, 
    lambda_hat, beta_hat, year_params, ATT, se, t_stat, p_val, Y, X, N, T, N0, T0, coef_table
  )
  return obj
end

function Base.summary(obj::SynthDID)
  
  Y_col = obj.Y_col
  D_col = obj.D_col
  se_method = obj.se_method
  coef_table = obj.coef_table
  covariates = obj.covariates
  cov_method = obj.cov_method

  println("Outcome variable: $Y_col")
  println("Treatment variable: $D_col")
  println("Standard error estimation method: $se_method")

  if !isnothing(covariates)
    print("Controled for covariate(s) ")
    for i in covariates
      print("$i, ")
    end
    print("using $cov_method residuals method\n")
  end

  println(coef_table)

end

function plot(obj::SynthDID, plottype::String = "outcomes")

  if !(plottype in ["outcomes", "weights"])
    throw(ArgumentError("positional argument plottype must be either \"outcomes\" or \"weights\""))
  end

  if plottype == "weights"
    res = Dict(
      "tyears" => obj.tyears, 
      "weights" => Dict("omega" => obj.omega_hat, "lambda" => obj.lambda_hat)
    )
    return plot_weights(res)
  end

  res = Dict(
    "t_span" => obj.t_span, "tyears" => obj.tyears, 
    "weights" => Dict("omega" => obj.omega_hat, "lambda" => obj.lambda_hat), 
    "year_params" => obj.year_params, "Y" => obj.Y
  )

  return plot_outcomes(res)
end

function Base.show(io::IO, ::MIME"text/html", obj::SynthDID)

  table = obj.coef_table
  println(io, table)
end
  
function Base.show(io::IO, ::MIME"text/plain", obj::SynthDID)

  table = obj.coef_table
  println(io, table)
end