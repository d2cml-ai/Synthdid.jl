function sparsify_function(v::Vector)
  v[v .<= maximum(v) / 4] .= 0
  return v ./ sum(v)
end

function data_setup(
  data::DataFrame, S_col::Union{String, Symbol}, 
  T_col::Union{String, Symbol}, D_col::Union{String, Symbol}
)

  tdf = copy(data)
  select!(groupby(tdf, S_col), :, D_col => maximum => :tunit)
  tdf.ty = @. ifelse(tdf[:, D_col] == 0, nothing, tdf[:, T_col])
  select!(groupby(tdf, S_col), :, :ty => minnothing => :tyear)
  sort!(tdf, [T_col, :tunit, S_col])

  return tdf
end

function projected(data, Y_col, S_col, T_col, covariates)

  k = size(covariates, 1)
  X = Matrix(data[:, covariates])
  y = data[:, Y_col]

  # Pick non-treated
  df_c = data[isnothing.(data.tyear), :]

  # One-hot encoding for T_col and S_col
  select!(df_c, :, [S_col => ByRow(isequal(v)) => Symbol(v) for v in unique(df_c[:, S_col])[2:end]])
  select!(df_c, :, [T_col => ByRow(isequal(v)) => Symbol(v) for v in unique(df_c[:, T_col])[2:end]])
  o_h_cov = Symbol.([covariates; unique(df_c[:, S_col])[2:end]; unique(df_c[:, T_col])[2:end]])

  # Create X_c Matrix with covariates, one-hot encoding for T_col and S_col. Create Y_c vector
  y_c = df_c[:, Y_col]
  X_c = Matrix(df_c[:, o_h_cov])

  # OLS of Y_c on X_c, get β
  XX = [X_c ones(size(X_c, 1))]' * [X_c ones(size(X_c, 1))]
  Xy = [X_c ones(size(X_c, 1))]' * y_c
  all_beta = inv(XX) * Xy
  beta = all_beta[1:k]

  # Calculate adjusted Y
  Y_adj = y - X * beta
  
  # Output projected data
  data[:, Y_col] = Y_adj
  return data, beta, X
end

function minnothing(x)
  x = x[.!isnothing.(x)]
  if length(x) == 0
    return nothing
  end
  return minimum(x)
end

function find_treat(W)
  N1 = 0
  for row in eachrow(W)
    if 1 in row
      N1 += 1
    end
  end
  return N1
end

function collapse_form(Y::Matrix, N0::Int64, T0::Int64)
  N, T = size(Y)
  Y_T0N0 = Y[1:N0, 1:T0]
  Y_T1N0 = mean(Y[1:N0, T0 + 1:end], dims = 2)
  Y_T0N1 = mean(Y[N0 + 1:end, 1:T0], dims = 1)
  Y_T1N1 = mean(Y[N0 + 1:end, T0 + 1:end])

  return [Y_T0N0 Y_T1N0; Y_T0N1 Y_T1N1]
end

# function pairwise_sum_decreasing(x::Vector{Number}, y::Vector{Number})
function pairwise_sum_decreasing(x::Vector, y::Vector)
  na_x = isnan.(x)
  na_y = isnan.(y)
  x[na_x] .= minimum(x[.!na_x])
  y[na_y] .= minimum(y[.!na_y])
  pairwise_sum = x .+ y
  pairwise_sum[na_x.&na_y] .= NaN
  return pairwise_sum
end

mutable struct random_walk
  Y::Matrix
  n0::Number
  t0::Number
  L::Matrix
end

function random_low_rank()
  n0 = 100
  n1 = 10
  t0 = 120
  t1 = 20
  n = n0 + n1
  t = t0 + t1
  tau = 1
  sigma = 0.5
  rank = 2
  rho = 0.7
  var = [rho^(abs(x - y)) for x in 1:t, y in 1:t]
  W = Int.(1:n .> n0) * transpose(Int.(1:t .> t0))

  # U = rand(Poisson(sqrt.(1:n) ./ sqrt(n)), n, rank)
  pU = Poisson(sqrt(sample(1:n)) ./ sqrt(n))
  pV = Poisson(sqrt(sample(1:t)) ./ sqrt(t))
  U = rand(pU, n, rank)
  V = rand(pV, t, rank)

  # sample.(1:n)

  alpha = reshape(repeat(10 * (1:n) ./ n, outer=(t, 1)), n, t)
  beta = reshape(repeat(10 * (1:t) ./ t, outer=(n, 1)), n, t)
  mu = U * V' + alpha + beta
  error = rand(pV, size(mu))
  Y = mu .+ tau .* W .+ sigma .* error
  random_data = random_walk(Y, n0, t0, mu)
  return random_data
end