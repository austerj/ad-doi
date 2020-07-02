struct Strangle <: AbstractContract
    T::Float64
    K₁::Float64
    K₂::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::Strangle)
    @unpack T, K₁, K₂ = contract
    τ = T-t

    put = u(t, s, σ, r, BlackScholes(), EuropeanPut(T,K₁))
    call = u(t, s, σ, r, BlackScholes(), EuropeanCall(T,K₂))

    put + call
end

function h(s, contract::Strangle)
    @unpack K₁,K₂ = contract
    max(K₁-s, 0) + max(s-K₂, 0)
end
