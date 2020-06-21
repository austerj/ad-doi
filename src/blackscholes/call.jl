struct EuropeanCall <: AbstractContract
    T::Float64
    K::Float64
end

@muladd function u(state::BlackScholesState, contract::EuropeanCall)
    @unpack t, s, σ, r = state
    @unpack T, K = contract
    τ = T-t

    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    u = exp(r*τ)*s*Φ(d1) - K*Φ(d2)
end

function h(s, contract::EuropeanCall)
    @unpack K = contract
    h = max(s-K, 0)
end
