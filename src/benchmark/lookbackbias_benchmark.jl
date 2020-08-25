include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.

nsteps_set = [2^k for k = 4:21]
nsteps_set_doi = nsteps_set[1:10]
nestimates = 1

contract = FloatingLookbackPut(T);

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

mc_price, doi_price = mc_benchmark(nsteps_set, 10^5, heston, contract, nestimates)

invsqrt_mc = invsqrtols(nsteps_set,mc_price[1][:])
ols_doi = log2linols(nsteps_set,doi_price[1][:])

mc_plot = plot(nsteps_set, invsqrt_mc, xaxis=:log2, alpha=.8, label="monte carlo", color=colors[1])
plot!(nsteps_set, ols_doi, xaxis=:log2, alpha=.8, label="doi", color=colors[2])
plot!(nsteps_set, doi_price[1], xaxis=:log2, alpha=.8, label="", seriestype=:scatter, color=colors[2])
plot!(nsteps_set, mc_price[1], xaxis=:log2, alpha=.8, label="", seriestype=:scatter, color=colors[1])

doi_plot = plot(nsteps_set, ols_doi, xaxis=:log2, alpha=.8, label="doi", color=colors[2])
plot!(nsteps_set, doi_price[1], xaxis=:log2, alpha=.8, label="", seriestype=:scatter, color=colors[2])

plot(mc_plot, doi_plot, layout=(1,2))
xlabel!("number of steps")
ylabel!("price estimate")
savefig("LookbackStepsize.pdf")
