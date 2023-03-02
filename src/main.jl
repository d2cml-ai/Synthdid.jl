mutable struct synthdid_est1
  estimate::Float64
  weight::Any
  setup::Any
  opts::Any
  N0::Int64
  T0::Int64
end

function sparsify_function(v::Vector)
  v[v.<=maximum(v)/4] .= 0
  return v ./ sum(v)
end

function synthdid_estimate(Y::Matrix, N0::Int, T0::Int;
  X::Array=Array{Any,3}(undef, size(Y, 1), size(Y, 2), 0),
  # X = zeros([size(Y, 1), size(Y, 2), 0])
  noise_level::Float64=std(diff(Y[1:N0, 1:T0], dims=2)),
  eta_omega::Float64=((size(Y, 1) - N0) * (size(Y, 2) - T0))^(1 / 4),
  eta_lambda::Float64=1e-6,
  zeta_omega::Float64=eta_omega * noise_level,
  zeta_lambda::Float64=eta_lambda * noise_level,
  omega_intercept::Bool=true,
  lambda_intercept::Bool=true,
  weights=Dict("omega" => nothing, "lambda" => nothing, "vals" => [1, 2, 3.0]),
  update_omega::Bool=isnothing(weights["omega"]),
  update_lambda::Bool=isnothing(weights["lambda"]),
  min_decrease::Float64=1e-5 * noise_level,
  max_iter::Int=10000,
  sparsify::Function=sparsify_function,
  max_iter_pre_sparsify::Int=100)
  if (!(size(Y)[1] > N0) && !(size(Y)[2] > T0) && !(length(size(X)) == 2 || length(size(X)) == 3) & !(size(X)[1:2] == size(Y)) && !(isa(weights, Dict))
      && !((isnothing(weights["lambda"])) || (length(weights["lambda"]) == T0)) && !((isnothing(weights["omega"])) || (length(weights["omega"]) == N0))
      && !(!(isnothing(weights["lambda"])) || (update_lambda)) && !((!isnothing(weights["omega"]) || (update_omega))))

    error("error at !(size(Y)[1] > N0) || !(size(Y)[2] > T0) || ... in synthdid_estimate function")
  else
    "continue"
  end

  N1 = size(Y, 1) - N0
  T1 = size(Y, 2) - T0

  if length(size(X)) == 3
    weights["vals"] = nothing
    weights["lambda_vals"] = nothing
    weights["omega_vals"] = nothing
    if (update_lambda)
      Yc = collapse_form(Y, N0, T0)
      lambda_opt = sc_weight_fw(
        # TODO: Lambda weights make a error
        Yc[1:N0, :], zeta_lambda, intercept=lambda_intercept, lambda=weights["lambda"], min_decrease=min_decrease, max_iter=max_iter_pre_sparsify)
      if (!isnothing(sparsify))
        lambda_lambda_opt = sparsify(lambda_opt["lambda"])
        lambda_opt = sc_weight_fw(Yc[1:N0, :], zeta_lambda, intercept=lambda_intercept, lambda=lambda_lambda_opt, min_decrease=min_decrease, max_iter=max_iter)
      end
      weights["lambda"] = lambda_opt["lambda"]
      weights["lambda_vals"] = lambda_opt["vals"]
      weights["vals"] = lambda_opt["vals"]
    end
    if (update_omega)
      Yc = collapse_form(Y, N0, T0)
      matrix_yc = Matrix(transpose(Yc[:, 1:T0]))
      omega_opt = sc_weight_fw(
        matrix_yc, zeta_omega, intercept=omega_intercept, lambda=weights["omega"], min_decrease=min_decrease, max_iter=max_iter_pre_sparsify)
      if (!isnothing(sparsify))
        omega_lambda_opt = sparsify(omega_opt["lambda"])
        omega_opt = sc_weight_fw(matrix_yc, zeta_omega, intercept=omega_intercept, lambda=omega_lambda_opt, min_decrease=min_decrease, max_iter=max_iter)
      end
      weights["omega"] = omega_opt["lambda"]
      weights["omega_vals"] = omega_opt["vals"]
      if (isnothing(weights["vals"]))
        weights["vals"] = omega_opt["vals"]
      else
        weights["vals"] = pairwise_sum_decreasing(weights["vals"], omega_opt["vals"])
      end
    end
  else
    YC = collapse_form(Y, N0, T0)
    Xc = [collapse_form(Xi, N0, T0) for Xi in eachslice(X, dims=3)]
    Yc = collapsed.form(Y, N0, T0)
  end

  weights["beta"] = nothing
  X_beta = contract3(X, weights["beta"])
  a = vcat(-weights["omega"], fill(1 / N1, N1))'
  b = vcat(-weights["lambda"], fill(1 / T1, T1))

  estimate = a * (Y .- X_beta) * b
  setup = Dict("Y" => Y, "X" => X, "N0" => N0, "T0" => T0)
  opts = Dict(
    "zeta_omega" => zeta_omega,
    "zeta_lambda" => zeta_lambda,
    "omega_intercept" => omega_intercept,
    "lambda_intercept" => lambda_intercept,
    "update_omega" => update_omega,
    "update_lambda" => update_lambda,
    "min_decrease" => min_decrease,
    "max_iter" => max_iter
  )
  return synthdid_est1(estimate, weights, setup, opts, N0, T0)
end




function sc_estimate(Y, N0, T0, eta_omega=1e-6; kargs...)

  estimate = synthdid_estimate(Y, N0, T0, eta_omega=1e-16, omega_intercept=false,
    weights=Dict("omega" => nothing, "lambda" => fill(0, T0), "vals" => [1, 2, 3.0]))
  return estimate
end



function did_estimate(Y, N0, T0; kargs...)
  estimate = synthdid_estimate(Y, N0, T0, weights=Dict("omega" => fill(1 / N0, N0), "lambda" => fill(1 / T0, T0), "vals" => [1, 2, 3.0]), kargs...)
  return estimate
end




# TODO: synthdid_placebo
function synthdid_placebo(estimate::synthdid_est1, terated_fraction=nothing)
  setup = estimate.setup
  opts = estimate.opts
  weights = estimate.weight
  x_beta = contract3(setup["X"], weights["beta"])
  estimator = estimate.estimate

  if (isnothing(terated_fraction))
    terated_fraction = 1 - setup["T0"] / size(setup.Y, 2)
  end
  placebo_t0 = floor(setup["T0"] * (1 - terated_fraction))
end

function synthdid_effect_curve(estimate::synthdid_est1)
  setup = estimate.setup
  weights = estimate.weight
  x_beta = contract3(setup["X"], weights["beta"])

  N1 = size(setup["Y"], 1) - estimate.N0
  T1 = size(setup["Y"], 2) - estimate.T0

  tau_sc = vcat(-weights["omega"], fill(1 / N1, N1))' * (setup["Y"] .- x_beta)
  tau_curve = tau_sc[setup["T0"].+(1:T1)] .- (tau_sc[1:setup["T0"]]' * weights["lambda"])
  return tau_curve
end

