struct SmoothedPowerBinaryPut <: AbstractContract
    T::Float64
    K::Float64
    n::Float64
    ε::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::SmoothedPowerBinaryPut)
    @unpack T, K, n, ε = contract
    τ = T-t

    d = (log(s/K) + (r+(n-0.5)*σ^2)*τ) / (√τ*σ)

    # smoothing term
    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    factor = K^n / ε
    
    dε1 = (log(s/(K+ε)) + (r+0.5*σ^2)*τ) / (√τ*σ)
    dε2 = dε1 - σ*√τ
    gap_term = K*Φ(-dε2) - exp(r*τ)*s*Φ(-dε1)
    put_term = K*Φ(-d2) - exp(r*τ)*s*Φ(-d1)

    smoothing_term = factor*(put_term - gap_term)

    exp((r+(n-1)*0.5*σ^2)*τ*n)*s^n*Φ(-d) + smoothing_term
end

function h(s, contract::SmoothedPowerBinaryPut)
    @unpack K, n = contract
    # does not include smoothing so that mc estimates can be benchmark for bias
    (s<K) * s^n
end
