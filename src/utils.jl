


function collapse_form(Y::Union{DataFrame,Matrix}, N0::Int64, T0::Int64)
  N, T = size(Matrix(Y))
  head = Y[1:N0, 1:T0]
  head_row_mean = mean(Y[1:N0, (T0+1):T], dims=2)
  head_matrix = hcat(head, head_row_mean)
  bottom_col_mean = mean(Y[(N0+1):N, 1:T0], dims=1)
  bottom = mean(Y[(N0+1):N, (T0+1):T])
  bottom_matrix = hcat(bottom_col_mean, bottom)
  return vcat(head_matrix, bottom_matrix)
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

mutable struct panelMatrix
  Y::Array
  W::Array
  names::Vector
  time::Vector
  N0::Int
  T0::Int
end

function panel_matrices(panel::DataFrame;
  unit::Union{String,Int}=1, time::Union{String,Int}=2,
  outcome::Union{String,Int}=3, treatment::Union{String,Int}=4,
  treated_last::Bool=true)
  index_to_name(x) = x in 1:size(panel, 2) ? names(panel)[x] : x
  if any(ismissing.(eachrow(panel)))
    error("Missing values in `panel`.")
  end
  keep = [unit, time, outcome, treatment]
  panel = panel[:, keep]
  panel = sort(panel, [unit, time])

  unique_years = unique(panel[:, time])
  unique_units = unique(panel[:, unit])

  num_years = length(unique(panel[:, time]))
  num_units = length(unique(panel[:, unit]))

  Y = reshape(panel[:, outcome], num_years, num_units)'
  W = reshape(panel[:, treatment], num_years, num_units)'

  w = sum(W, dims=2)
  w = [x > 0 for x in w]
  N0 = num_units - sum(w)
  treat_time = any(W, dims=1)
  T0 = ([i for i in 1:length(treat_time) if treat_time[i]] |> first) - 1

  unit_order = if treated_last
    sortperm(W[:, T0+1])
  else
    collect(1:size(Y, 1))
  end
  panel = panelMatrix(Y[unit_order, :], W[unit_order, :], unique_units, unique_years, N0, T0)
  return panel
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