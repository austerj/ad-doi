include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.
contract = EuropeanCall(T,K);

nsteps = 200
npaths_set = [Int(floor(10^k)) for k = 1:0.25:6]
nestimates = 50
heston_price = 15.9400786350804

mc, doi = mc_benchmark(nsteps, npaths_set, heston, contract, nestimates, heston_price)

theme(
    :default,
    size=(800, 300),
    xgrid=false,
    ygrid=false,
    markerstrokewidth=0,
    markersize=5,
    linewidth=1.5
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]

mc_fit = log10ols(npaths_set, mc[1])
doi_fit = log10ols(npaths_set, doi[1])

plot(npaths_set, mc_fit, xaxis=:log10, yaxis=:log10, alpha=.8, label="monte carlo", color=colors[1])
plot!(npaths_set, mc[1], xaxis=:log10, yaxis=:log10, alpha=.8, label="", seriestype=:scatter, color=colors[1])
plot!(npaths_set, doi[1], xaxis=:log10, yaxis=:log10, alpha=.8, label="", seriestype=:scatter, color=colors[2])
plot!(npaths_set, doi_fit, xaxis=:log10, yaxis=:log10, alpha=.8, label="doi", color=colors[2])
xlabel!("number of paths")
ylabel!("mean relative error")

savefig("MonteCarloBenchmark.pdf")
