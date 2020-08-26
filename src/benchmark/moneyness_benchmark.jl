include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Plots
using ColorSchemes
using DelimitedFiles

pyplot()
heston_prices = readdlm("heston_prices.csv", ',', Float64)

# contract
T = 1.
Ks = 20:5:180.

nsteps = 200
npaths = 10
nestimates = 10000

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
    ygrid=false,
    markerstrokewidth=0,
    markersize=5,
    linewidth=1.5
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]

mc_fit = log10linols(Ks, mc_re[1,:])
doi_fit = log10linols(log10.(Ks), doi_re[1,:])

plot(Ks, mc_fit, yaxis=:log10, alpha=.8, label="monte carlo", color=colors[1])
plot!(Ks, mc_re[1,:], yaxis=:log10, alpha=.8, label="", seriestype=:scatter, color=colors[1])
plot!(Ks, doi_fit, yaxis=:log10, alpha=.8, label="doi", color=colors[2])
plot!(Ks, doi_re[1,:], yaxis=:log10, alpha=.8, label="", seriestype=:scatter, color=colors[2])
xlabel!("K")
ylabel!("mean relative error")
savefig("MoneynessRelativeErrors.pdf")
