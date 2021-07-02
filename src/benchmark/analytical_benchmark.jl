include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes

pyplot()

# contract
T = 1.
K = 100.

@muladd function analytical_diffop(t, s, ν, model::Heston, contract::EuropeanCall, state::AbstractState)::Float64
    @unpack r, κ, θ, ξ, ρ = model
    @unpack T, K = contract

    τ = T-t

    if τ == 0
        return zero(s)
    end

    # σ̄ and derivatives wrt ν
    σ̄ = √(θ + (ν-θ)/κ*(1-exp(-κ*τ))/τ)
    σ̄∂ν = 0.5/σ̄*(1-exp(-κ*τ))/(κ*τ)
    σ̄∂νν = -σ̄∂ν^2/σ̄

    d1 = (log(s/K) + (r+0.5*σ̄^2)*τ) / (√τ*σ̄)
    d2 = d1 - σ̄*√τ

    # sensitivities
    ∂νν = exp(r*τ)*s*φ(d1)*(d1*d2-1)/σ̄^3*(1-exp(-κ*τ))^2/(4*κ^2*τ^(3/2))
    ∂sν = -exp(r*τ)*φ(d1)*d2/σ̄^2*(1-exp(-κ*τ))/(2*κ*τ)

    ξ*ν*(ρ*s*∂sν + 0.5*ξ*∂νν)
end

function analytical_estimator(nsteps, npaths, model::Heston, contract::AbstractContract, rng::AbstractRNG)
    @unpack s₀, ν₀, r = model
    @unpack T = contract

    Δ = T/nsteps
    t = Δ:Δ:T
    initial_state = state(model, contract)

    u₀ = u(0., s₀, ν₀, model, contract, initial_state)
    ad_A₀ = Δ*diffop(0., s₀, ν₀, model, contract, initial_state)/2
    an_A₀ = Δ*analytical_diffop(0., s₀, ν₀, model, contract, initial_state)/2

    nthreads = Threads.nthreads()
    npath_atomic = Threads.Atomic{Int}(0)
    payoff_atomic = Threads.Atomic{Float64}(0)

    payoffs = Vector{Float64}(undef, npaths)
    ad_doi = Vector{Float64}(undef, npaths)
    an_doi = Vector{Float64}(undef, npaths)

    for j = 1:npaths
        s, ν = s₀, ν₀

        ad_A = ad_A₀
        an_A = an_A₀
        path_state = deepcopy(initial_state)

        for i = 1:nsteps
            s, ν = step(s, ν, Δ, model, rng)
            state!(t, s, ν, contract, path_state)
            ad_A += diffop(t[i], s, ν, model, contract, path_state)
            an_A += analytical_diffop(t[i], s, ν, model, contract, path_state)
        end
        # trapezoidal rule; diffop is zero at time T hence division not needed
        ad_A *= Δ
        an_A *= Δ

        payoffs[j] = h(s, contract, path_state)
        ad_doi[j] = ad_A
        an_doi[j] = an_A
    end

    discount = exp(-r*T)
    payoffs *= discount
    ad_doi = discount * (u₀ .+ ad_doi)
    an_doi = discount * (u₀ .+ an_doi)

    payoffs, ad_doi, an_doi
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
colors = [ColorSchemes.ice[i] for i in 64:64:192]

# contract
T = 1.
K = 100
contract = EuropeanCall(T,K)

nsteps = 200
npaths = 100

rng = Xorshift128Plus(1)
payoffs, ad_doi, an_doi = analytical_estimator(nsteps, npaths, heston, contract, rng)

payoff_ad = plot(ad_doi, payoffs, alpha=.8, seriestype=:scatter, color=colors[1], xlabel="analytical doi", ylabel="payoff", legend=false)
payoff_an = plot(an_doi, payoffs, alpha=.8, seriestype=:scatter, color=colors[2], xlabel="dual doi", ylabel="payoff", legend=false)
ad_an = plot(ad_doi, an_doi, alpha=.8, seriestype=:scatter, color=colors[3], xlabel="dual doi", ylabel="analytical doi", legend=false)
plot(payoff_ad, payoff_an, ad_an, layout=(1,3))
savefig("AnalyticalDOI.pdf")
