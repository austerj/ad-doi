include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Plots
using ColorSchemes
using DelimitedFiles

pyplot()
heston_prices = readdlm("heston_prices.csv", ',', Float64)
# mc_mean = readdlm("../mc_price_means.csv", ',', Float64)
# mc_std = readdlm("../mc_price_stds.csv", ',', Float64)
# doi_mean = readdlm("../doi_price_means.csv", ',', Float64)
# doi_std = readdlm("../doi_price_stds.csv", ',', Float64)

# contract
T = 1.
Ks = 60:5:140.

nsteps = 50
npaths = 10
nestimates = 1000

# price benchmark
mc = Array{Float64, 2}(undef, 2, length(Ks))
doi = Array{Float64, 2}(undef, 2, length(Ks))

for (j,K) in enumerate(Ks)
    contract = EuropeanCall(T,K);
    mc_tuple, doi_tuple = mc_benchmark(nsteps, npaths, heston, contract, nestimates)
    mc[1,j] = mc_tuple[1][1]
    mc[2,j] = mc_tuple[2][1]
    doi[1,j] = doi_tuple[1][1]
    doi[2,j] = doi_tuple[2][1]
end

# relative errors benchmark
mc_re = Array{Float64, 2}(undef, 2, length(Ks))
doi_re = Array{Float64, 2}(undef, 2, length(Ks))

for (j,K) in enumerate(Ks)
    contract = EuropeanCall(T,K);
    mc_tuple, doi_tuple = mc_benchmark(nsteps, npaths, heston, contract, nestimates, heston_prices[j])
    mc_re[1,j] = mc_tuple[1][1]
    mc_re[2,j] = mc_tuple[2][1]
    doi_re[1,j] = doi_tuple[1][1]
    doi_re[2,j] = doi_tuple[2][1]
end

theme(
    :default,
    size=(800, 300),
    xgrid=false,
    ygrid=false
)
colors = [ColorSchemes.ice[i] for i in 64:64:192]

mc_plot = plot(Ks, mc[1,:], ribbon=mc[2,:], label="monte carlo", alpha=.8, fillapha=.5, color=colors[1], xlabel="K")
heston_plot = plot(Ks, heston_prices, ribbon=zeros(length(Ks)), label="numerical integration", alpha=.8, fillapha=.5, color=colors[2], xlabel="K")
doi_plot = plot(Ks, doi[1,:], ribbon=doi[2,:], label="doi", alpha=.8, fillapha=.5, color=colors[3], xlabel="K")

plot(heston_plot, doi_plot, mc_plot, layout=(1,3))
savefig("MoneynessPrices.pdf")

mc_re_plot = plot(Ks, mc_re[1,:], ribbon=mc_re[2,:], label="monte carlo", alpha=.8, fillapha=.5, color=colors[1], xlabel="K")
doi_re_plot = plot(Ks, doi_re[1,:], ribbon=doi_re[2,:], label="doi", alpha=.8, fillapha=.5, color=colors[3], xlabel="K")

plot(doi_re_plot, mc_re_plot, layout=(1,2))
savefig("MoneynessRelativeErrors.pdf")
