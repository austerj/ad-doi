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
