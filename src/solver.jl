
# When fw_step and sc_weight_fw are used with omega, the columns must represent units and rows must represent time.
# The opposite is true for lambda
function fw_step(
  A::Matrix, b::Vector, x::Vector; eta::Number, 
  alpha::Union{Nothing,Float64} = nothing
)::Vector{Float64}
  Ax = A * x
  half_grad = (Ax .- b)' * A + eta * x'
  i = findmin(half_grad)[2][2]
  if !isnothing(alpha)
    x *= (1 - alpha)
    x[i] += alpha
    return x
  else
    d_x = -x
    d_x[i] = 1 - x[i]
    if all(d_x .== 0)
      return x
    end
    d_err = A[:, i] - Ax
    step_upper = -half_grad * d_x
    step_bot = sum(d_err .^ 2) + eta * sum(d_x .^ 2)
    step = step_upper[1] / step_bot
    constrained_step = min(1, max(0, step))
    return x + constrained_step * d_x
  end
end

function sc_weight_fw(
  A::Matrix, b::Vector, x::Union{Vector, Nothing} = nothing; 
  intercept::Bool = true, zeta::Number,
  min_decrease::Number = 1e-3, max_iter::Int64 = 1000
)
  
  k = size(A, 2)
  n = size(A, 1)
  if isnothing(x)
    x = fill(1 / k, k)
  end
  if intercept
    A = A .- mean(A, dims = 1)
    b = b .- mean(b, dims = 1)
  end

  t = 0
  vals = zeros(max_iter)
  eta = n * real(zeta ^ 2)
  while (t < max_iter) && (t < 2 || vals[t-1] - vals[t] > min_decrease ^ 2)
    t += 1
    x_p = fw_step(A, b, x, eta = eta)
    x = x_p
    err = A * x - b
    vals[t] = real(zeta^2) * sum(x .^ 2) + sum(err .^ 2) / n
  end
  Dict("params" => x, "vals" => vals)
end

# sc_weight_covariates only runs when there are covariates specified
# This implements the procedure as in Abadie et al. (2010)
# X::Vector{Matrix{Float64}} w/ covariates
# Y (outcome) and X must come from unstacked data (year columns, unit rows, covariate values) and must be in collapsed form
# For staggered treatments, this function should be applied for a in A (see Clarke et al. 2023, p. 9-10)

function sc_weight_covariates(
  Y::Matrix, X::Vector; zeta_lambda = 0, zeta_omega = 0, 
  lambda_intercept::Bool = true, omega_intercept::Bool = true, 
  min_decrease::Float64 = 1e-3, max_iter::Int = 1000, 
  lambda = nothing, omega = nothing, beta = nothing, 
  update_lambda::Bool = true, update_omega::Bool = true
)

  T0 = size(Y, 2) - 1
  N0 = size(Y, 1) - 1

  if isnothing(lambda) lambda = fill(1 / T0, T0) end
  if isnothing(omega) omega = fill(1 / N0, N0) end
  if isnothing(beta) beta = zeros(size(X, 1)) end

  function update_weights(Y, lambda, omega)
    
    Y_lambda = if lambda_intercept Y[1:N0, :] .- mean(Y[1:N0, :], dims = 1) else Y[1:N0, :] end
    if update_lambda lambda = fw_step(Y_lambda[:, 1:T0], Y_lambda[:, T0 + 1], lambda, eta = N0 * real(zeta_lambda ^ 2)) end
    err_lambda = Y_lambda * [lambda; -1]

    Y_omega = if omega_intercept Y'[1:T0, :] .- mean(Y'[1:T0, :], dims = 1) else Y[:, 1:T0]' end
    if update_omega omega = fw_step(Y_omega[:, 1:N0], Y_omega[:, N0 + 1], omega, eta = T0 * real(zeta_omega ^ 2)) end
    err_omega = Y_omega * [omega; -1]

    val = real(zeta_omega ^ 2) * sum(omega .^ 2) + real(zeta_lambda ^ 2) * sum(lambda .^ 2) + sum(err_omega .^ 2) / T0 + sum(err_lambda .^ 2) / N0

    return Dict("val" => val, "lambda" => lambda, "omega" => omega, "err_lambda" => err_lambda, "err_omega" => err_omega)
  end

  vals = zeros(max_iter)
  t = 0
  Y_beta = Y - sum(beta .* X)
  weights = update_weights(Y_beta, lambda, omega)

  while (t < max_iter) && ((t < 2) || (vals[t - 1] - vals[t] > min_decrease ^ 2))
    
    t += 1
    gr_lambda = (Ref(weights["err_lambda"]') .* [arr[1:N0, :] for arr in X]) .* Ref([weights["lambda"]; -1]) ./ N0
    gr_omega = (Ref(weights["err_omega"]') .* [arr[:, 1:T0]' for arr in X]) .* Ref([weights["omega"]; -1]) ./ T0
    grad_beta = -(gr_lambda[1] + gr_omega[1])

    alpha = 1 / t
    beta = beta .- alpha * grad_beta
    Y_beta = Y - sum(beta .* X)
    weights = update_weights(Y_beta, weights["lambda"], weights["omega"])
    vals[t] = weights["val"]
  end
  return Dict("lambda" => weights["lambda"], "omega" => weights["omega"], "beta" => beta, "vals" => vals)
end


