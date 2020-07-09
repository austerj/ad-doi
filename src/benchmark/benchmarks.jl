function npaths_benchmark(nsteps, npaths, heston, contract, nestimates)
    mc_estimates = Vector{Float64}(undef, nestimates)
    doi_estimates = Vector{Float64}(undef, nestimates)

    for j = 1:nestimates
        mc, doi = parallel_estimator(nsteps, npaths, heston, contract)
        mc_estimates[j] = mc
        doi_estimates[j] = doi
    end

    mc_estimates, doi_estimates
end

function mc_benchmark(nsteps_set, npaths_set, heston, contract, nestimates)
    mc_mean = Array{Float64,2}(undef, length(nsteps_set), length(npaths_set))
    mc_std = Array{Float64,2}(undef, length(nsteps_set), length(npaths_set))
    doi_mean = Array{Float64,2}(undef, length(nsteps_set), length(npaths_set))
    doi_std = Array{Float64,2}(undef, length(nsteps_set), length(npaths_set))

    for (j, npaths) in enumerate(npaths_set)
        for (i, nsteps) in enumerate(nsteps_set)
            npaths = npaths_set[j]
            mc_estimates, doi_estimates = npaths_benchmark(nsteps, npaths, heston, contract, nestimates)

            mc_mean[i,j] = mean(mc_estimates)
            mc_std[i,j] = std(mc_estimates)
            doi_mean[i,j] = mean(doi_estimates)
            doi_std[i,j] = std(doi_estimates)
        end
    end
    
    (mc_mean, mc_std), (doi_mean, doi_std)
end

function mc_benchmark(nsteps, npaths_set, heston, contract, nestimates, heston_price)
    mc_re_mean = Vector{Float64}(undef, length(npaths_set))
    mc_re_std = Vector{Float64}(undef, length(npaths_set))
    doi_re_mean = Vector{Float64}(undef, length(npaths_set))
    doi_re_std = Vector{Float64}(undef, length(npaths_set))

    for (j, npaths) in enumerate(npaths_set)
        npaths = npaths_set[j]
        mc_estimates, doi_estimates = npaths_benchmark(nsteps, npaths, heston, contract, nestimates)

        mc_re = abs.((mc_estimates .- heston_price)) / heston_price
        doi_re = abs.((doi_estimates .- heston_price)) / heston_price

        mc_re_mean[j] = mean(mc_re)
        mc_re_std[j] = std(mc_re)
        doi_re_mean[j] = mean(doi_re)
        doi_re_std[j] = std(doi_re)
    end
    
    (mc_re_mean, mc_re_std), (doi_re_mean, doi_re_std)
end
