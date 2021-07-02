include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.
contract = EuropeanCall(T,K);

@noinline function dual_sensitivities(t, s, ν, model::Heston, contract::AbstractContract)
    # add small pertubations to prevent caching
    s += randn()*1e-5
    ν += randn()*1e-5

    # create dual numbers
    zs = Dual(Dual(s,1.,0.), Dual(0.,0.,0.))
    zν = Dual(Dual(ν,0.,1.), Dual(1.,0.,0.))

    # σ̄ computed inside u
    uz = u(t, zs, zν, model, contract)

    ∂s = uz.value.partials[1]  # delta
    ∂ν = uz.value.partials[2]  # vega
    ∂sν = uz.partials[1].partials[1]  # vanna
    ∂νν = uz.partials[1].partials[2]  # vomma

    # force return to prevent compiler from skipping computations
    ∂s, ∂ν, ∂sν, ∂νν
end

@noinline @muladd function analytical_sensitivities(t, s, ν, model::Heston, contract::EuropeanCall)
    @unpack r, κ, θ = model
    @unpack T, K = contract

    # add small pertubations to prevent caching
    s += randn()*1e-5
    ν += randn()*1e-5

    # σ̄ and derivatives wrt ν
    τ = T-t
    σ̄ = √(θ + (ν-θ)/κ*(1-exp(-κ*τ))/τ)
    σ̄∂ν = 0.5/σ̄*(1-exp(-κ*τ))/(κ*τ)
    σ̄∂νν = -σ̄∂ν^2/σ̄

    # sensitivities
    d1 = (log(s/K) + (r+0.5*σ̄^2)*τ) / (√τ*σ̄)
    d2 = d1 - σ̄*√τ
    ∂s = exp(r*τ)*Φ(d1)  # delta
    ∂ν = exp(r*τ)*s*φ(d1)*√τ*σ̄∂ν  # vega
    ∂νν = exp(r*τ)*s*φ(d1)*(d1*d2-1)/σ̄^3*(1-exp(-κ*τ))^2/(4*κ^2*τ^(3/2))  # vomma
    ∂sν = -exp(r*τ)*φ(d1)*d2/σ̄^2*(1-exp(-κ*τ))/(2*κ*τ)  # vanna

    # force return to prevent compiler from skipping computations
    ∂s, ∂ν, ∂sν, ∂νν
end

nsamples = 10000
nevals = 1000
nseconds = 120

T = 1
K = 100.
contract = EuropeanCall(T,K);

dual_benchmark = @benchmark dual_sensitivities(0., $s₀, $ν₀, $heston, $contract) samples=nsamples seconds=nseconds evals=nevals
analytical_benchmark = @benchmark analytical_sensitivities(0., $s₀, $ν₀, $heston, $contract) samples=nsamples seconds=nseconds evals=nevals

dt = dual_benchmark.times
at = analytical_benchmark.times

theme(
    :default,
    size=(800, 300),
    palette=[ColorSchemes.ice[i] for i in 64:128:192],
    xgrid=false,
    xminorticks=2
)

nbins = 100
thresh = .75
alpha = 0.8

range_start = max(min(mean(dt), mean(at)) - thresh*max(std(dt), std(at)), 0)
range_stop = max(mean(dt), mean(at)) + thresh*max(std(dt), std(at))
timerange = range(range_start, stop=range_stop, length=nbins)

histogram(at, bins=timerange, normalize=:probability, fillalpha=alpha, label="analytical")
histogram!(dt, bins=timerange, normalize=:probability, fillalpha=alpha, label="dual")
xlabel!("nanoseconds")
ylabel!("relative frequency")
yaxis!(bordercolor="white")

savefig("DualBenchmark.pdf")
