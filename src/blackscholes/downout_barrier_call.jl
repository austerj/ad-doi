using SpecialFunctions: zeta

struct DownOutBarrierCall <: AbstractContract
    T::Float64
    K::Float64
    H::Float64
    nsteps::Int32
end

# track breaching of barrier
mutable struct BarrierBreach <: AbstractState
    barrier_breached::Bool
end

function state(model::AbstractModel, contract::DownOutBarrierCall)
    BarrierBreach(model.s₀ <= contract.H)
end

function state!(t, s, ν, contract::AbstractContract, state::BarrierBreach)
    state.barrier_breached = state.barrier_breached || s <= contract.H
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::DownOutBarrierCall, state::BarrierBreach)
    @unpack T, K, H, nsteps = contract
    @unpack barrier_breached = state

    if barrier_breached || s < contract.H
        return zero(s)
    end

    # barrier bias correction for discretely-monitored price
    β = -zeta(.5) / √(2*π)
    H = H * exp(-β*σ*√(T/nsteps))

    τ = T-t

    # case dependent variables
    if H < K
        d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
        h1 = (log(H^2/(s*K)) + (r+0.5*σ^2)*τ) / (√τ*σ)
    else
        d1 = (log(s/H) + (r+0.5*σ^2)*τ) / (√τ*σ)
        h1 = (log(H/s) + (r+0.5*σ^2)*τ) / (√τ*σ)
    end

    d2 = d1 - σ*√τ
    h2 = h1 - σ*√τ
    p = 2*r/σ^2

    exp(r*τ)*s*(Φ(d1)-Φ(h1)*(H/s)^(p+1)) - K*(Φ(d2)-Φ(h2)*(H/s)^(p-1))
end

function h(s, contract::DownOutBarrierCall, state::BarrierBreach)
    state.barrier_breached ? zero(s) : max(s - contract.K, 0)
end

