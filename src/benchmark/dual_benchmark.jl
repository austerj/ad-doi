include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.
contract = EuropeanCall(T,K);

@noinline function dualsensitivities(s, ν, model::Heston, contract::AbstractContract)
    # add small pertubations to prevent caching
    s += randn()*1e-5
    ν += randn()*1e-5

    zs = Dual(Dual(s,1.,0.), Dual(0.,0.,0.))
    zν = Dual(Dual(ν,0.,1.), Dual(1.,0.,0.))

    # σ̄ computed inside u
    uz = u(0., zs, zν, model, contract)

    ∂s = uz.value.partials[1]  # delta
    ∂ν = uz.value.partials[2]  # vega
    ∂sν = uz.partials[1].partials[1]  # vanna
    ∂νν = uz.partials[1].partials[2]  # vomma

    # force return to prevent compiler from skipping computations
    ∂s, ∂ν, ∂sν, ∂νν
end

@noinline @muladd function analyticalsensitivities(s, ν, model::Heston, contract::AbstractContract)
    @unpack r, κ, θ = model
    @unpack T, K = contract

    # add small pertubations to prevent caching
    s += randn()*1e-5
    ν += randn()*1e-5

    # σ̄ and derivatives wrt ν
    σ̄ = √((θ + (ν-θ)/κ*(1-exp(-κ*T))/T))
    σ̄∂ν = 0.5*(1-exp(-κ*T))/√(κ*T*(κ*θ+(ν-θ)*(1-exp(-κ*T))))
    σ̄∂νν = -0.25*(1-exp(-(κ*T)))^2 / (√(κ*T)*(κ*θ+(ν-θ)*(1-exp(-κ*T)))^(3/2))

    d1 = (log(s/K) + (r+0.5*σ̄^2)*T) / (√T*σ̄)
    d2 = d1 - σ̄*√T

    # Black-Scholes sensitivities
    BS_∂ν = exp(r*T)*s*φ(d1)*√T  # vega
    BS_∂sν = -exp(r*T)*φ(d1)*(d2/σ̄)  # vanna
    BS_∂νν = exp(r*T)*s*φ(d1)*√T*(d1*d2)/σ̄  # vomma

    ∂s = exp(r*T)*Φ(d1)  # delta
    ∂ν = BS_∂ν*σ̄∂ν  # vega
    ∂sν = BS_∂sν*σ̄∂ν  # vanna
    ∂νν = BS_∂νν*σ̄∂ν^2 + BS_∂ν*σ̄∂νν  # vomma

    # force return to prevent compiler from skipping computations
    ∂s, ∂ν, ∂sν, ∂νν
end

nsamples = 10000
nevals = 1000
nseconds = 120

dualbenchmark = @benchmark dualsensitivities($s₀, $ν₀, $heston, $contract) samples=nsamples seconds=nseconds evals=nevals
analyticalbenchmark = @benchmark analyticalsensitivities($s₀, $ν₀, $heston, $contract) samples=nsamples seconds=nseconds evals=nevals

dt = dualbenchmark.times
at = analyticalbenchmark.times

theme(
    :default,
    size=(800, 300),
    palette=[ColorSchemes.ice[i] for i in 64:128:192],
    xgrid=false,
    xminorticks=2
)

nbins = 100
thresh = 0.5
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
