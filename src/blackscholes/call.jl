struct EuropeanCall <: AbstractContract
    T::Float64
    K::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::EuropeanCall)
    @unpack T, K = contract
    τ = T-t

    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    exp(r*τ)*s*Φ(d1) - K*Φ(d2)
end

function h(s, contract::EuropeanCall)
    @unpack K = contract
    max(s-K, 0)
end
