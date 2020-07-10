include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

rng = Xoroshiro128Plus(1)

# contract
T = 1.
K = 100.
n = 2.

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

function mc_benchmark(nsteps, npaths_set, heston, contract, nestimates)
    mc_mean = Vector{Float64}(undef, length(npaths_set))
    mc_std = Vector{Float64}(undef, length(npaths_set))
    doi_mean = Vector{Float64}(undef, length(npaths_set))
    doi_std = Vector{Float64}(undef, length(npaths_set))

    for (j, npaths) in enumerate(npaths_set)
        npaths = npaths_set[j]
        mc_estimates, doi_estimates = npaths_benchmark(nsteps, npaths, heston, contract, nestimates)

        mc_mean[j] = mean(mc_estimates)
        mc_std[j] = std(mc_estimates)
        doi_mean[j] = mean(doi_estimates)
        doi_std[j] = std(doi_estimates)
    end
    
    (mc_mean, mc_std), (doi_mean, doi_std)
end

nsteps = 200
npaths_set = [Int(floor(10^k)) for k = 1:0.25:5]
nestimates = 50

contract = PowerBinaryPut(T, K, n);
mc_price, doi_price = mc_benchmark(nsteps, npaths_set, heston, contract, nestimates)

theme(
    :default,
    size=(800, 300),
    palette=[ColorSchemes.ice[i] for i in 64:128:192],
    xgrid=false,
    ygrid=false
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]

plot(npaths_set, mc_price[1], ribbon=mc_price[2], yaxis=:log10, xaxis=:log10, label="monte carlo", alpha=.8, fillapha=.5)
plot!(npaths_set, doi_price[1], ribbon=doi_price[2], yaxis=:log10, xaxis=:log10, label="doi", alpha=.8, fillapha=.5)
plot(npaths_set, mc_price[1], ribbon=mc_price[2], xaxis=:log10, label="monte carlo", alpha=.8, fillapha=.5)
plot!(npaths_set, doi_price[1], ribbon=doi_price[2], xaxis=:log10, label="doi", alpha=.8, fillapha=.5)
xlabel!("number of paths")
ylabel!("price estimate")

savefig("PowerBinaryPutBenchmark.pdf")
