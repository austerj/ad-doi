struct EuropeanPut <: AbstractContract
    T::Float64
    K::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::EuropeanPut)
    @unpack T, K = contract
    τ = T-t

    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    K*Φ(-d2) - exp(r*τ)*s*Φ(-d1)
end

function h(s, contract::EuropeanPut)
    @unpack K = contract
    max(K-s, 0)
end
