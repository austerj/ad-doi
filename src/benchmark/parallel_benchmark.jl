include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

rng = Xoroshiro128Plus(1)

# contract
T = 1.
K = 100.
contract = EuropeanCall(T,K);

nsteps = 50
npaths = 200

nsamples = 2000
nevals = 5
nseconds = 120

function create_heston()
    s₀ = 100
    ν₀ = 0.16
    r = 0.04
    κ = 0.6
    θ = 0.04
    ξ = 0.2
    ρ = -0.15

    # add small pertubations to prevent caching
    s₀ += randn()*1e-5
    ν₀ += randn()*1e-5
    model = Heston(s₀, ν₀, r, κ, θ, ξ, ρ)
end

function estimator_benchmark(nsteps, npaths, contract, rng)
    model = create_heston()
    mc, doi = estimator(nsteps, npaths, model, contract, rng)
end

function parallel_estimator_benchmark(nsteps, npaths, contract)
    model = create_heston()
    mc, doi = parallel_estimator(nsteps, npaths, model, contract)
end

singlebenchmark = @benchmark estimator_benchmark($nsteps, $npaths, $contract, $rng) samples=nsamples seconds=nseconds evals=nevals
parallelbenchmark = @benchmark parallel_estimator_benchmark($nsteps, $npaths, $contract) samples=nsamples seconds=nseconds evals=nevals

# convert to microseconds
st = singlebenchmark.times ./ 1000
pt = parallelbenchmark.times ./ 1000

theme(
    :default,
    size=(800, 300),
    palette=[ColorSchemes.ice[i] for i in 64:128:192],
    xgrid=false,
    xminorticks=2
)

nbins = 120
thresh = 3.5
alpha = 0.8

range_start = max(min(mean(pt), mean(st)) - thresh*max(std(pt), std(st)), 0)
range_stop = max(mean(pt), mean(st)) + thresh*max(std(pt), std(st))
timerange = range(range_start, stop=range_stop, length=nbins)

histogram(st, bins=timerange, normalize=:probability, fillalpha=alpha, label="single-threaded")
histogram!(pt, bins=timerange, normalize=:probability, fillalpha=alpha, label="multi-threaded")
xlabel!("microseconds")
ylabel!("relative frequency")
yaxis!(bordercolor="white")

savefig("ParallelBenchmark.pdf")
