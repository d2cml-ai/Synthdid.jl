

function ar2_correlation_matrix(ar_coef, T)
    result = repeat([0.0], T)
    result[1] = 1
    result[2] = ar_coef[1] / (1 - ar_coef[2])

    for t in 3:T
        result[t] = ar_coef[1] * result[t-1] + ar_coef[2] * result[t-2]
    end

    index_matrix = abs.((1:T)' .- (1:T)'') .+ 1
    cor_matrix = reshape(result[index_matrix], (T, T))

    return cor_matrix
end

function fit_ar2(E)
    T_full = size(E)[2]
    E_ts = E[:, 3:T_full]
    E_lag_1 = E[:, 2:(T_full-1)]
    E_lag_2 = E[:, 1:(T_full-2)]

    a_1 = sum(diag(E_lag_1 * E_lag_1'))
    a_2 = sum(diag(E_lag_2 * E_lag_2'))
    a_3 = sum(diag(E_lag_1 * E_lag_2'))

    matrix_factor = hcat([a_1, a_3], [a_3, a_2])
    b_1 = sum(diag(E_lag_1 * E_ts'))
    b_2 = sum(diag(E_lag_2 * E_ts'))

    ar_coef = inv(matrix_factor) * [b_1, b_2]

    return ar_coef
end

function randomize_treatment(pi, N, N1)
    if length(pi) == 1
        p = Binomial.(1, pi)
        assignment_sim = rand(p, N)

    elseif N <= length(pi)
        pi1 = pi[1:N]
        assignment_sim = []
        for i in collect(1:N)
            p = Binomial.(1, pi1[i])
            assignment_sim = vcat(assignment_sim, rand(p, 1))
        end
    elseif N > length(pi)
        pi1 = repeat(vec(pi), N ÷ length(pi) + 1)
        assignment_sim = []
        for i in collect(1:N)
            p = Binomial.(1, pi1[i])
            assignment_sim = vcat(assignment_sim, rand(p, 1))
        end
    end
    index_as = findall(assignment_sim .== 1)
    if sum(assignment_sim) > N1
        index_pert = sample(index_as, N1)
        assignment_sim = repeat([0.0], N)
        assignment_sim[index_pert] .= 1
    elseif sum(assignment_sim) == 0
        index_pert = sample(1:N, N1)
        assignment_sim = repeat([0.0], N)
        assignment_sim[index_pert] .= 1
    end

    return assignment_sim
end

function lindsey_density_estimate(x, K, deg)
    x_min = minimum(x)
    x_max = maximum(x)

    range_x = x_max - x_min
    low_x = x_min - 0.2 * range_x                             # 20% step outside of range of x
    up_x = x_max + 0.2 * range_x
    range_full = up_x - low_x
    splits = collect(range(low_x, up_x, length=K + 1))       # split the range into K segments
    mesh_size = splits[2] - splits[1]
    centers = (splits[Not(1)] .+ splits[Not(K + 1)]) ./ 2
    counts = vec(freqtable(cut(x, splits)))
    scale = sum(counts) * mesh_size

    data_matrix1 = ns(centers, df=deg)
    data_matrix = DataFrame(data_matrix1, :auto)
    data_matrix[:, "counts"] = counts
    function form(n, y="counts", x_v="x")
        y1 = StatsModels.Term(Symbol(y))
        eq = []
        eq = vcat(eq, StatsModels.Term(Symbol("x1")))
        vl = 2
        while vl < n + 1
            eq = vcat(eq, eq[vl-1] + StatsModels.Term(Symbol(x_v * string(vl))))
            vl = vl + 1
        end
        return y1 ~ eq[n]
    end
    form2 = form(size(data_matrix1)[2])
    pois_reg_res = glm(form2, data_matrix, Poisson())
    counts_pois = exp.(hcat([1 1 1 1 1]', data_matrix1) * coef(pois_reg_res))
    dens_pois = counts_pois ./ scale

    return Dict("centers" => centers, "density" => dens_pois)
end

function decompose_Y(Y, j_rank)
    N = size(Y)[1]
    T = size(Y)[2]

    svd_data_mat = svd(Y)
    factor_unit = Matrix(svd_data_mat.U[:, 1:j_rank] * sqrt(N))
    factor_time = Matrix(svd_data_mat.V[:, 1:j_rank] * sqrt(T))

    magnitude = svd_data_mat.S[1:j_rank] ./ sqrt(N * T)
    L = factor_unit * Matrix(Diagonal(magnitude)) * factor_time'

    E = Y - L
    F = mean(L, dims=2) .* repeat([1.0], T)' .+ repeat([1.0], N)[:, :] * mean(L, dims=1) .- mean(L)
    M = L .- F

    return Dict("F" => F, "M" => M, "E" => E, "unit_factors" => factor_unit)
end

function estimate_dgp(Y, assignment_vector, j_rank)
    N = size(Y)[1]
    T = size(Y)[2]
    overall_mean = mean(Y)
    overall_sd = norm(Y .- overall_mean) / sqrt(N * T)
    Y_norm = (Y .- overall_mean) ./ overall_sd

    components = decompose_Y(Y_norm, j_rank)
    M = components["M"]
    F = components["F"]
    E = components["E"]
    unit_factors = components["unit_factors"]

    ar_coef = round.(fit_ar2(E), digits=2)
    cor_matrix = ar2_correlation_matrix(ar_coef, T)
    scale_sd = norm(E' * E ./ N) ./ norm(cor_matrix)
    cov_mat = cor_matrix .* scale_sd

    data = DataFrame(unit_factors[:, :], :auto)
    data[:, "assignment_vector"] = vec(assignment_vector)
    function form(n, y="assignment_vector", x_v="x")
        y1 = StatsModels.Term(Symbol(y))
        eq = []
        eq = vcat(eq, StatsModels.Term(Symbol("x1")))
        vl = 2
        while vl < n + 1
            eq = vcat(eq, eq[vl-1] + StatsModels.Term(Symbol(x_v * string(vl))))
            vl = vl + 1
        end
        return y1 ~ eq[n]
    end
    form2 = form(size(data)[2] - 1)
    assign_prob = predict(glm(form2, data, Binomial()))

    return Dict("F" => F, "M" => M, "Sigma" => cov_mat, "pi" => assign_prob, "ar_coef" => ar_coef)
end

function simulate_dgp(parameters, N1, T1)
    F = parameters["F"]
    M = parameters["M"]
    Sigma = parameters["Sigma"]
    pi = parameters["pi"]

    N = size(M)[1]
    T = size(M)[2]

    assignment = randomize_treatment(pi, N, N1)
    # assignment = [1 1 0 0]
    N1 = sum(assignment)
    N0 = N - N1
    T0 = T - T1

    Y = F + M + rand(MvNormal(Sigma), N)'
    data = DataFrame(assignment=vec(assignment), id=1:length(vec(assignment)))
    sort_id = sort(data)[:, "id"]
    return Dict("Y" => Y[sort_id], "N0" => N0, "T0" => T0)
end

