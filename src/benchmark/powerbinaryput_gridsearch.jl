include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.
n = 2.

nsteps_set = [2^k for k = 4:10]
npaths_set = [Int(floor(10^k)) for k = 1:0.5:6]
nestimates = 50

contract = PowerBinaryPut(T, K, n);
mc_price, doi_price = mc_benchmark(nsteps_set, npaths_set, heston, contract, nestimates)

theme(
    :default,
    size=(800, 300),
    # palette=[ColorSchemes.ice[i] for i in 64:128:192],
    # color=:ice,
    xgrid=false,
    ygrid=false,
    colorbar=true
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]

## load csv files
# using DelimitedFiles
# mc_mean = readdlm("../mc_price_means.csv", ',', Float64)
# mc_std = readdlm("../mc_price_stds.csv", ',', Float64)
# doi_mean = readdlm("../doi_price_means.csv", ',', Float64)
# doi_std = readdlm("../doi_price_stds.csv", ',', Float64)

mc_price = (mc_mean, mc_std)
doi_price = (doi_mean, doi_std)

# plot mean estimates from 10^6 samples
mc_plot = plot(nsteps_set, mc_price[1][:,end], ribbon=mc_price[2][:,end],
    xaxis=:log2, label="monte carlo", alpha=.8, fillapha=.5, color=colors[1],
    xlabel="number of steps", ylabel="price estimate")
doi_plot = plot(nsteps_set, doi_price[1][:,end], ribbon=doi_price[2][:,end],
    xaxis=:log2, label="doi", alpha=.8, fillapha=.5, color=colors[2],
    xlabel="number of steps", ylabel="price estimate")
plot(doi_plot, mc_plot)
savefig("PowerBinaryPutStepsize.pdf")
