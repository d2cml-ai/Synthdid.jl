
function sdid(
  data, Y_col::Union{String, Symbol}, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol};
  covariates::Union{Vector{String}, Vector{Symbol}, Nothing} = nothing, 
  cov_method = "optimized", noise_level::Union{Float64, Nothing} = nothing, 
  eta_omega::Union{Float64, Nothing} = nothing, eta_lambda::Union{Float64, Nothing} = 1e-6,
  zeta_omega::Union{Float64, Nothing} = nothing, zeta_lambda::Union{Float64, Nothing} = nothing, 
  omega_intercept::Bool = true, lambda_intercept::Bool = true,
  min_decrease::Union{Float64, Nothing} = nothing, max_iter::Int = 10000,
  sparsify::Union{Function, Nothing} = sparsify_function, 
  max_iter_pre_sparsify::Int = 100, vce = "placebo"
)

  att = Float64[]

  Y_col = Symbol(Y_col)
  S_col = Symbol(S_col)
  T_col = Symbol(T_col)
  D_col = Symbol(D_col)

  if all(in.(["tunit", "ty", "tyear"], Ref(names(data))))
    tdf = copy(data)
  else
    tdf = data_setup(data, S_col, T_col, D_col)
  end

  tyears = sort(unique(tdf.tyear)[.!isnothing.(unique(tdf.tyear))])
  T_total = sum(Matrix(unstack(tdf, S_col, T_col, D_col)[:, 2:end]))
  units = unique(data[:, S_col])
  N_out = size(units, 1)

  tdf_ori = copy(tdf)
  info_names = ["year", "tau", "weighted_tau", "N0", "T0", "N1", "T1"]
  year_params = DataFrame([[] for i in info_names], info_names)
  T_out = size(unique(tdf[:, T_col]), 1)
  info_beta = []

  # project covariates method
  if !isnothing(covariates) && cov_method == "projected"
    covariates = Symbol.(covariates)
    for covariate in covariates
      year_params[:, "beta_" * String(covariate)] = []
    end

    # create DataFrame with projected data
    tdf, info_beta = projected(tdf, Y_col, S_col, T_col, covariates)
  end

  # estimation when no covariates or already projected covariates
  if isnothing(covariates) || cov_method == "projected"

    # for a in A
    for year in tyears
      info = []
      df_y = tdf[in.(tdf.tyear, Ref([year, nothing])), [Y_col, S_col, T_col, :tunit]]
      N1 = size(unique(df_y[df_y.tunit .== 1, S_col]), 1)
      T1 = maximum(data[:, T_col]) - year + 1
      T_post = N1 * T1
      # create Y matrix and collapse it
      Y = Matrix(unstack(df_y, S_col, T_col, Y_col)[:, 2:end])
      N, T = size(Y)
      N0 = N - N1
      T0 = T - T1
      Yc = collapse_form(Y, N0, T0)

      # calculate penalty parameters
      noise_level = std(diff(Y[1:N0, 1:T0], dims = 2)) # gotta fix this, probably its own function
      eta_omega = ((size(Y, 1) - N0) * (size(Y, 2) - T0))^(1 / 4)
      eta_lambda = 1e-6
      zeta_omega = eta_omega * noise_level
      zeta_lambda = eta_lambda * noise_level
      min_decrease = 1e-5 * noise_level

      # calculate lambda
      lambda_opt = sc_weight_fw(
        Yc[1:N0, 1:T0], Yc[1:N0, end], nothing, 
        intercept = lambda_intercept, 
        zeta = zeta_lambda, 
        min_decrease = min_decrease, 
        max_iter = max_iter_pre_sparsify
      )
      
      if !isnothing(sparsify)
        lambda_opt = sc_weight_fw(
          Yc[1:N0, 1:T0], Yc[1:N0, end], sparsify(lambda_opt["params"]), 
          intercept = lambda_intercept, 
          zeta = zeta_lambda, 
          min_decrease = min_decrease, 
          max_iter = max_iter
        )
      end

      lambda = lambda_opt["params"]

      # calculate omega
      omega_opt = sc_weight_fw(
        Yc'[1:T0, 1:N0], Yc[end, 1:T0], nothing, 
        intercept = omega_intercept, 
        zeta = zeta_omega, 
        min_decrease = min_decrease, 
        max_iter = max_iter_pre_sparsify
      )
      
      if !isnothing(sparsify)
        omega_opt = sc_weight_fw(
          Yc'[1:T0, 1:N0], Yc[end, 1:T0], sparsify(omega_opt["params"]), 
          intercept = omega_intercept, 
          zeta = zeta_omega, 
          min_decrease = min_decrease, 
          max_iter = max_iter
        )
      end

      omega = omega_opt["params"]

      # calculate tau for this a
      tau_hat = [-omega; fill(1/N1, N1)]' * Y * [-lambda; fill(1/T1, T1)]
      tau_w = T_post / T_total * tau_hat
      att = [att; tau_w]
      info = [year tau_hat tau_w N0 T0 N1 T1]
      if !isnothing(covariates)
        info = [info info_beta']
      end
      info_df = DataFrame([i for i in info], names(year_params))
      append!(year_params, info_df)
    end
    
    att = sum(att)
  end

  # optimized errors with covariates method
  if !isnothing(covariates) && cov_method == "optimized"
    covariates = Symbol.(covariates)
    for covariate in covariates
      year_params[:, "beta_" * String(covariate)] = []
    end

    # for a in A
    for year in tyears
      info = []
      df_y = tdf[in.(tdf.tyear, Ref([year, nothing])), [[Y_col, S_col, T_col, :tunit]; covariates]]
      N1 = size(unique(df_y[df_y.tunit .== 1, S_col]), 1)
      T1 = maximum(data[:, T_col]) - year + 1
      T_post = N1 * T1

      # create Y matrix
      Y = Matrix(unstack(df_y, S_col, T_col, Y_col)[:, 2:end])
      N = size(Y, 1)
      T = size(Y, 2)
      N0 = N - N1
      T0 = T - T1
      Yc = collapse_form(Y, N0, T0)

        # create penalty parameters
      noise_level = std(diff(Y[1:N0, 1:T0], dims = 2)) # this needs fixed, maybe in its own function
      eta_omega = ((size(Y, 1) - N0) * (size(Y, 2) - T0))^(1 / 4)
      eta_lambda = 1e-6
      zeta_omega = eta_omega * noise_level
      zeta_lambda = eta_lambda * noise_level
      min_decrease = 1e-5 * noise_level

      # create X vector of matrices
      X = []
      for covar in covariates
        X_temp = Matrix(unstack(df_y, S_col, T_col, covar)[:, 2:end])
        push!(X, X_temp)
      end

      Xc = collapse_form.(X, N0, T0)
      
      # find lambda, omega, and beta for Yc with X
      weights = sc_weight_covariates(
        Yc, Xc, zeta_lambda = zeta_lambda, zeta_omega = zeta_omega, 
        lambda_intercept = lambda_intercept, omega_intercept = omega_intercept, 
        min_decrease = min_decrease, max_iter = max_iter, lambda = nothing, 
        omega = nothing
      )

      # calculate tau for this a
      X_beta = sum(weights["beta"] .* X)
      tau_hat = [-weights["omega"]; fill(1/N1, N1)]' * (Y - X_beta) * [-weights["lambda"]; fill(1/T1, T1)]
      tau_w = T_post / T_total * tau_hat
      att = [att; tau_w]
      info = [year tau_hat tau_w N0 T0 N1 T1]
      info = [info weights["beta"]']
      info_df = DataFrame([i for i in info], names(year_params))
      append!(year_params, info_df)
    end
    
    # weighed average of tau
    att = sum(att)
  end
  
  out = Dict("att" => att, "year_params" => year_params, "T" => T_out, "N" => N_out, "data" => data, "proc_data" => tdf_ori)

  return out
end

# function sc_estimate(Y, N0, T0, eta_omega=1e-6; kargs...)

#   estimate = synthdid_estimate(Y, N0, T0, eta_omega=1e-16, omega_intercept=false,
#     weights=Dict("omega" => nothing, "lambda" => fill(0, T0), "vals" => [1, 2, 3.0]))
#   return estimate
# end



# function did_estimate(Y, N0, T0; kargs...)
#   estimate = synthdid_estimate(Y, N0, T0, weights=Dict("omega" => fill(1 / N0, N0), "lambda" => fill(1 / T0, T0), "vals" => [1, 2, 3.0]), kargs...)
#   return estimate
# end




# TODO: synthdid_placebo
# function synthdid_placebo(estimate::synthdid_est1, terated_fraction=nothing)
#   setup = estimate.setup
#   opts = estimate.opts
#   weights = estimate.weight
#   x_beta = contract3(setup["X"], weights["beta"])
#   estimator = estimate.estimate

#   if (isnothing(terated_fraction))
#     terated_fraction = 1 - setup["T0"] / size(setup.Y, 2)
#   end
#   placebo_t0 = floor(setup["T0"] * (1 - terated_fraction))
# end

# function synthdid_effect_curve(estimate::synthdid_est1)
#   setup = estimate.setup
#   weights = estimate.weight
#   x_beta = contract3(setup["X"], weights["beta"])

#   N1 = size(setup["Y"], 1) - estimate.N0
#   T1 = size(setup["Y"], 2) - estimate.T0

#   tau_sc = vcat(-weights["omega"], fill(1 / N1, N1))' * (setup["Y"] .- x_beta)
#   tau_curve = tau_sc[setup["T0"].+(1:T1)] .- (tau_sc[1:setup["T0"]]' * weights["lambda"])
#   return tau_curve
# end

