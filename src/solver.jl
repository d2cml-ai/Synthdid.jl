
function contract3(X::Array{Float64,3}, v::Union{Vector,Nothing}=nothing)::Matrix{Float64}
  if !isnothing(v) && size(X, 3) != length(v)
    throw(ArgumentError("The length of `v` must match the size of the third dimension of `X`"))
  end
  out = zeros(Float64, size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end
function contract3(X, v::Union{Vector,Nothing}=nothing)::Matrix{Float64}
  if !isnothing(v) && size(X, 3) != length(v)
    throw(ArgumentError("The length of `v` must match the size of the third dimension of `X`"))
  end
  out = zeros(Float64, size(X, 1), size(X, 2))
  if isnothing(v)
    return out
  end
  for ii in eachindex(v)
    out .+= v[ii] * X[:, :, ii]
  end
  return out
end

function fw_step(A::Matrix, x::Vector{Float64}; b::Vector{Float64}, eta::Number, alpha::Union{Nothing,Float64}=nothing)::Vector{Float64}
  Ax = A * x
  half_grad = transpose(Ax .- b) * A + (eta * x)'
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
  Y::Matrix, zeta::Number;
  intercept::Bool=true, lambda::Union{Vector,Nothing}=nothing,
  min_decrease::Number=1e-3, max_iter::Int64=1000)
  T0 = size(Y, 2) - 1
  N0 = size(Y, 1)
  if isnothing(lambda)
    lambda = fill(1 / T0, T0)
  end
  if intercept
    Y = Y .- mean(Y, dims=1)
  end

  t = 0
  vals = zeros(max_iter)
  A = Y[:, 1:T0]
  b = Y[:, T0+1]
  eta = N0 * real(zeta^2)
  while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
    t += 1
    lambda_p = fw_step(A, lambda, b=b, eta=eta)
    lambda = lambda_p
    err = Y[1:N0, :] * [lambda; -1]
    vals[t] = real(zeta^2) * sum(lambda .^ 2) + sum(err .^ 2) / N0
  end
  Dict("lambda" => lambda, "vals" => vals)
end;

mutable struct sc_weight_fw_covariates1
  lambda
  omega
  beta
  vals
end

mutable struct update_weights1
  val
  lambda
  omega
  err_lambda
  err_omega
end

function sc_weight_fw_covariates(Y::Matrix; X=cat(zeros(size(Y)), dims=ndims(Y) + 1),
  zeta_lambda=0, zeta_omega=0,
  lambda_intercept=true, omega_intercept=true,
  min_decrease=1e-3, max_iter=1000,
  lambda=nothing, omega=nothing, beta=nothing, update_lambda=true, update_omega=true)

  if length(size(Y)) == 2 && length(size(X)) == 3 && all(size(Y) == size(X)[1:2]) && all(Matrix((isfinite.(Y)))) && all(isfinite.(X))
    "continue"
  else
    error("the following condition is not met: length(size(Y)) != 2 || length(size(X)) != 3 || all(size(Y) != size(X)[1:2]) || !all(Matrix((isfinite.(Y)))) || !all(isfinite.(X))")
  end

  T0 = size(Y)[2] - 1
  N0 = size(Y)[1] - 1
  if ndims(X) == 2
    cat(X; dims=ndims(X) + 1)
  end
  if isnothing(lambda)
    lambda = repeat([1 / T0], T0)
  end
  if isnothing(omega)
    omega = repeat([1 / N0], N0)
  end
  if isnothing(beta)
    beta = repeat([0.0], size(X)[3] - 1)
    if isempty(beta)
      beta = nothing
    end
  end

  function update_weights(Y, lambda, omega)

    Y_lambda = zeros(N0, T0 + 1)
    if lambda_intercept
      for i in 1:size(Y, 2)
        Y_lambda[:, i] = Y[1:N0, i] .- mean(Y[1:N0, i])
      end
    else
      Y_lambda = Y[1:N0, :]
    end
    if update_lambda
      lambda = fw_step(Y_lambda[:, 1:T0], lambda, b=Y_lambda[:, T0+1], eta=N0 * real(zeta_lambda^2))
    end
    err_lambda = Y_lambda * vcat(lambda, -1)

    Y_omega = zeros(size(Matrix(Y[:, 1:T0])', 1), size(Matrix(Y[:, 1:T0])', 2))
    if omega_intercept
      for i in 1:size(Matrix(Y[:, 1:T0])', 2)
        Y_omega[:, i] = Matrix(Y[:, 1:T0])'[1:T0, i] .- mean(Matrix(Y[:, 1:T0])'[1:T0, i])
      end
    else
      Y_omega = Matrix(Y[:, 1:T0])'
    end
    if update_omega
      omega = fw_step(Y_omega[:, 1:N0], omega, b=Y_omega[:, N0+1], eta=N0 * real(zeta_omega^2))
    end
    err_omega = Y_omega * vcat(omega, -1)
    val = real(zeta_omega .^ 2) * sum(omega .^ 2) + real(zeta_lambda .^ 2) * sum(lambda .^ 2) + sum(err_omega .^ 2) / T0 .+ sum(err_lambda .^ 2) ./ N0
    # return Dict("val" => val, "lambda" => lambda, "omega" => omega, "err_lambda" => err_lambda, "err_omega" => err_omega)
    res1 = update_weights1(val, lambda, omega, err_lambda, err_omega)
    return res1

  end

  vals = repeat([0.0], max_iter)
  t = 0
  Y_beta = Y .- contract3(X, beta)
  weights = update_weights(Y_beta, lambda, omega)

  while t < max_iter && (t < 2 || vals[t-1] - vals[t] > min_decrease^2)
    t = t + 1
    if size(X)[3] - 1 == 0
      grad_beta = 0
    else
      grad_beta = print("error in  while t < max_iter && (t < 2 || vals[t - 1] - vals[t] > min_decrease^2)")
    end

    alpha = 1 / t

    if isnothing(beta)
      beta = 0
    end

    beta = beta .- alpha * grad_beta
    beta = nothing
    Y_beta = Y .- contract3(X, beta)
    weights = update_weights(Y_beta, weights.lambda, weights.omega)
    vals[t] = weights.val
  end
  res2 = sc_weight_fw_covariates1(weights.lambda, weights.omega, beta, vals)

  return res2
end
