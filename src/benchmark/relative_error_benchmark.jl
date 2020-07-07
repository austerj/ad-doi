include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

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

# contract
T = 1.
K = 100.
contract = EuropeanCall(T,K);

nsteps = 50
npaths_set = [Int(floor(10^k)) for k = 1:0.25:6]
nestimates = 50
heston_price = 15.9400786350804

mc, doi = mc_benchmark(nsteps, npaths_set, heston, contract, nestimates, heston_price)

theme(
    :default,
    size=(800, 300),
    palette=[ColorSchemes.ice[i] for i in 64:128:192],
    xgrid=false,
    ygrid=false
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]

plot(npaths_set, mc[1], ribbon=mc[2], xaxis=:log10, yaxis=:log10, label="monte carlo", alpha=.8, fillapha=.5)
plot!(npaths_set, doi[1], ribbon=doi[2], xaxis=:log10, yaxis=:log10, label="doi", alpha=.8, fillapha=.5)
xlabel!("number of paths")
ylabel!("relative error")

savefig("MonteCarloBenchmark.pdf")
