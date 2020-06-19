struct EuropeanPut <: AbstractContract
    r::Float64
    T::Float64
    K::Float64
end

function u(state::BlackScholesState, contract::EuropeanPut)
    @unpack s, σ = state
    @unpack r, T, K = contract

    d1 = (log(s/K) + (r+0.5*σ^2)*T) / (√T*σ);
    d2 = d1 - σ*√T

    u = K*Φ(-d2) - exp(r*T)*s*Φ(-d1)
end

function h(state::BlackScholesState, contract::EuropeanPut)
    @unpack s = state
    @unpack K = contract

    h = max(K-s, 0)
end
