struct ExchangeOption <: AbstractContract
    T::Float64
end

@muladd function u(t, s, ν, r, ρ, model::BlackScholes, contract::ExchangeOption)
    @unpack T = contract
    τ = T-t

    σ = √(ν[1]^2 + ν[2]^2 - 2*ν[1]*ν[2]*ρ[1,2])

    d1 = (log(s[1]/s[2]) + 0.5*σ^2*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    exp(r*τ)*(s[1]*Φ(d1) - s[2]*Φ(d2))
end

function h(s, contract::ExchangeOption)
    max(s[1]-s[2], 0)
end
