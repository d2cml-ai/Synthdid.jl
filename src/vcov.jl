

function vcov_synthdid_estimate(object; method="bootstrap", replications=200)
    if method == "bootstrap"
        se = bootstrap_se(object, replications)
    elseif method == "jackknife"
        se = jackknife_se(object)
    elseif method == "placebo"
        se = placebo_se(object, replications)
    end
    se^2
end

function bootstrap_se(estimate, replications)
    sqrt((replications - 1) / replications) * std(bootstrap_sample(estimate, replications))
end

function bootstrap_sample(estimate, replications)
    setup = estimate.setup
    opts = estimate.opts
    weights = estimate.weight
    if setup["N0"] == size(setup["Y"])[1] - 1
        return NaN
    end
    function theta(ind)
        if all(ind .<= setup["N0"]) || all(ind .> setup["N0"])
            NaN
        else
            weights1 = copy(weights)
            weights_boot = weights1
            weights_boot["omega"] = sum_normalize(weights["omega"][sort(ind[ind.<=setup["N0"]])])
            synthdid_estimate(setup["Y"][sort(ind), :], sum(ind .<= setup["N0"]), setup["T0"], X=setup["X"][sort(ind), :, :], weights=weights_boot,
                zeta_omega=opts["zeta_omega"], zeta_lambda=opts["zeta_lambda"], omega_intercept=opts["omega_intercept"],
                lambda_intercept=opts["lambda_intercept"], update_omega=opts["update_omega"], update_lambda=opts["update_lambda"],
                min_decrease=opts["min_decrease"], max_iter=opts["max_iter"])
        end
    end

    bootstrap_estimates = repeat([NaN], replications)
    count = 0
    while count < replications
        a = theta(sample(1:size(setup["Y"])[1], size(setup["Y"])[1], replace=true))
        if typeof(a) == synthdid_est1
            bootstrap_estimates[count+1] = a.estimate
            count = count + 1
        end
    end
    bootstrap_estimates
end

function jackknife_se(estimate; weights=estimate.weight)
    setup = estimate.setup
    opts = estimate.opts
    weights = estimate.weight
    if !isnothing(weights)
        opts["update_omega"] = opts["update_lambda"] = false
    end
    if setup["N0"] == size(setup["Y"])[1] - 1 || (!isnothing(weights) && sum(weights["omega"] != 0) == 1)
        return NaN
    end

    function jackknife(x)
        n = length(x)
        u = repeat([0.0], n)
        function theta(ind)
            weights1 = copy(weights)
            weights_jk = weights1
            if !isnothing(weights)
                weights_jk["omega"] = sum_normalize(weights["omega"][ind[ind.<=setup["N0"]]])
            end
            synthdid_estimate(setup["Y"][ind, :], sum(ind .<= setup["N0"]), setup["T0"], X=setup["X"][ind, :, :], weights=weights_jk,
                zeta_omega=opts["zeta_omega"], zeta_lambda=opts["zeta_lambda"], omega_intercept=opts["omega_intercept"],
                lambda_intercept=opts["lambda_intercept"], update_omega=opts["update_omega"], update_lambda=opts["update_lambda"],
                min_decrease=opts["min_decrease"], max_iter=opts["max_iter"])
        end
        for i in 1:n
            u[i] = theta(x[Not(i)]).estimate
        end
        jack_se = sqrt(((n - 1) / n) * (n - 1) * var(u))

        jack_se
    end
    jackknife(1:size(setup["Y"])[1])
end

function placebo_se(estimate, replications)
    setup = estimate.setup
    opts = estimate.opts
    weights = estimate.weight
    N1 = size(setup["Y"])[1] - setup["N0"]
    if setup["N0"] <= N1
        error("must have more controls than treated units to use the placebo se")
    end
    function theta(ind)
        N0 = length(ind) - N1
        weights1 = copy(weights)
        weights_boot = weights1
        weights_boot["omega"] = sum_normalize(weights["omega"][ind[1:N0]])
        synthdid_estimate(setup["Y"][ind, :], N0, setup["T0"], X=setup["X"][ind, :, :], weights=weights_boot,
            zeta_omega=opts["zeta_omega"], zeta_lambda=opts["zeta_lambda"], omega_intercept=opts["omega_intercept"],
            lambda_intercept=opts["lambda_intercept"], update_omega=opts["update_omega"], update_lambda=opts["update_lambda"],
            min_decrease=opts["min_decrease"], max_iter=opts["max_iter"])
    end

    a = 0
    for i in 1:replications
        if i == 1
            a = theta(shuffle(1:setup["N0"])).estimate
        else
            a = vcat(a, theta(shuffle(1:setup["N0"])).estimate)
        end
    end
    sqrt((replications - 1) / replications) * std(a)
end

function sum_normalize(x)
    if sum(x) != 0
        x / sum(x)
    else
        repeat([1 / length(x)], length(x))
    end
end