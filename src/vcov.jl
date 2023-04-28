function jackknife_se(
  data, Y_col, S_col, T_col, D_col; att::Union{Float64, Nothing} = nothing, 
  omega_hat::Union{DataFrame, Nothing} = nothing, lambda_hat::Union{Dict, Nothing} = nothing, 
  kwargs...
)
  att_jk = []
  Y_col = Symbol(Y_col)
  S_col = Symbol(S_col)
  T_col = Symbol(T_col)
  units = unique(data[:, S_col])
  T_total = sum(Matrix(unstack(data, S_col, T_col, D_col)[:, 2:end]))

  if all(in.(["tunit", "ty", "tyear"], Ref(names(data))))
    sort!(data, [T_col, :tunit, S_col])
  else
    data = data_setup(data, S_col, T_col, D_col)
  end

  if isnothing(att) || isnothing(omega_hat) || isnothing(lambda_hat)
    res = sdid(data, Y_col, S_col, T_col, D_col; kwargs...)
    att = res["att"]
    omega_hat = res["weights"]["omega"]
    lambda_hat = res["weights"]["lambda"]
    if any(res["year_params"][:, "N1"] .< 2)
      throw(ErrorException("Jackknife standard error needs at least two treated units for each treatment period"))
    end
  end

  tyears = parse.(Int, names(omega_hat)[2:end])
  
  function theta(unit)

    # create aux data and aux omega
    aux_att = []
    aux_data = data[data[:, S_col] .!= unit, :]
    aux_omega = omega_hat[omega_hat[:, S_col] .!= unit, :]

    # for a in A
    for year in tyears
      # yearly params
      aux_lambda = lambda_hat[string(year)]
      df_y = aux_data[in.(aux_data.tyear, Ref([year, nothing])), [Y_col, S_col, T_col, :tunit]]
      N1 = size(unique(df_y[df_y.tunit .== 1, S_col]), 1)
      if N1 < 1 throw(ErrorException("Jackknife standard error needs at least two treated units for each treatment period")) end
      T1 = maximum(aux_data[:, T_col]) - year + 1
      T_post = N1 * T1

      # create Y matrix
      Y = Matrix(unstack(df_y, S_col, T_col, Y_col)[:, 2:end])

      # calculate tau for this a
      tau_hat = [-aux_omega; fill(1/N1, N1)]' * Y * [-aux_lambda; fill(1/T1, T1)]
      tau_w = T_post / T_total * tau_hat
      aux_att = [aux_att; tau_w]
    end

    # find aux att
    # aux_att = aux_res["att"]
    return sum(aux_att)
  end

  for unit in units
    aux_att = theta(unit)
    att_jk = [att_jk; aux_att]
  end

  # find jk_se
  jk_se = (N - 1)/N * sum((att_jk .- att) .^ 2)

  return jk_se ^ (1/2)
end

function bootstrap_se(data, Y_col, S_col, T_col, D_col; n_boot = 50, kwargs...)

  att_bs = []
  Y_col = Symbol(Y_col)
  S_col = Symbol(S_col)
  T_col = Symbol(T_col)
  units = unique(data[:, S_col])
  N = length(units)
  if all(in.(["tunit", "ty", "tyear"], Ref(names(data))))
    tdf = copy(data)
  else
    tdf = data_setup(data, S_col, T_col, D_col)
  end
  print("Bootstrap replications (",  n_boot, "). This may take some time.\n")
  print("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5\n")

  function theta(df)
    sample = rand(1:N, N)
    data_aux = copy(df)
    sampled_data = reduce(vcat, [groupby(data_aux, S_col)[i] for i in sample])
    
    if all(sampled_data[:, D_col] .== 1) || all(sampled_data[:, D_col] .== 0)
      theta(df)
    end
    
    for time in groupby(sampled_data, T_col)
      time[:, :id] = collect(1: size(time, 1))
    end
    boot_res = sdid(sampled_data, Y_col, :id, T_col, D_col; kwargs...)

    return boot_res["att"]
  end

  i = 1
  while i <= n_boot
    print(".")
    if i % 50 == 0
      print("     ", i, "\n")
    end
    aux_att = theta(tdf)
    att_bs = [att_bs; aux_att]
    i += 1
  end

  print("\n")
  bs_se = 1/n_boot * sum((att_bs .- sum(att_bs / n_boot)) .^ 2)

  return bs_se ^ (1/2)
end

function placebo_se(data, Y_col, S_col, T_col, D_col; n_reps = 50, kwargs...)
  
  att_p = []
  Y_col = Symbol(Y_col)
  S_col = Symbol(S_col)
  T_col = Symbol(T_col)

  if all(in.(["tunit", "ty", "tyear"], Ref(names(data))))
    tdf = copy(data)
  else
    tdf = data_setup(data, S_col, T_col, D_col)
  end

  tr_years = tdf[(.!isnothing.(tdf.tyear)) .&& (tdf[:, T_col] .== tdf.tyear), T_col]
  N_tr = size(tr_years, 1)
  df_co = tdf[tdf.tunit .== 0, :]
  N_co = size(unique(df_co[:, S_col]), 1)
  N_aux = N_co - N_tr
  print("Placebo replications (",  n_reps, "). This may take some time.\n")
  print("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5\n")
  
  function theta(df)

    aux_data = copy(df)
    aux_data = groupby(aux_data, S_col)[shuffle(1:end)]
  
    for i in 1:N_tr
      aux_data[N_aux + i][:, :tyear] .= tr_years[i]
    end

    aux_data = combine(aux_data, :)

    aux_data[(.!isnothing.(aux_data.tyear)) .&& (aux_data[:, T_col] .>= aux_data.tyear), D_col] .= 1
    aux_data = aux_data[:, Not(:tunit)]
    select!(groupby(tdf, S_col), :, D_col => maximum => :tunit)

    res_p = sdid(aux_data, Y_col, S_col, T_col, D_col; kwargs...)
    aux_att = res_p["att"]
    return aux_att
  end

  t = 1
  while t <= n_reps
    print(".")
    if t % 50 == 0
      print("     ", t, "\n")
    end
    aux_att = theta(df_co)
    att_p = [att_p; aux_att]
    t += 1
  end

  print("\n")
  p_se = 1/n_reps * sum((att_p .- sum(att_p / n_reps)) .^ 2)
  return p_se ^ (1/2)
end

# function sum_normalize(x)
#     if sum(x) != 0
#         x / sum(x)
#     else
#         repeat([1 / length(x)], length(x))
#     end
# end
