struct PowerBinaryPut <: AbstractContract
    T::Float64
    K::Float64
    n::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::PowerBinaryPut)
    @unpack T, K, n = contract
    τ = T-t

    d = (log(s/K) + (r+(n-0.5)*σ^2)*τ) / (√τ*σ)

    u = exp((r+(n-1)*0.5*σ^2)*τ*n)*s^n*Φ(-d)
end

function h(s, contract::PowerBinaryPut)
    @unpack K, n = contract
    h = (s<K) * s^n
end
