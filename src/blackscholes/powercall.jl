struct PowerCall <: AbstractContract
    T::Float64
    K::Float64
    n::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::PowerCall)
    @unpack T, K, n = contract
    τ = T-t

    d1 = (log(s/K^(1/n)) + (r+(n-0.5)*σ^2)*τ) / (√τ*σ)
    d2 = d1 - n*σ*√τ

    u = exp((r+(n-1)*0.5*σ^2)*τ*n)*s^n*Φ(d1) - K*Φ(d2)
end

function h(s, contract::PowerCall)
    @unpack K, n = contract
    h = max(s^n-K, 0)
end
