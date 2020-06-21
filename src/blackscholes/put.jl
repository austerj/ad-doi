struct EuropeanPut <: AbstractContract
    T::Float64
    K::Float64
end

@muladd function u(state::BlackScholesState, contract::EuropeanPut)
    @unpack t, s, σ, r = state
    @unpack T, K = contract
    τ = T-t

    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ);
    d2 = d1 - σ*√τ

    u = K*Φ(-d2) - exp(r*τ)*s*Φ(-d1)
end

function h(state::BlackScholesState, contract::EuropeanPut)
    @unpack s = state
    @unpack K = contract

    h = max(K-s, 0)
end
