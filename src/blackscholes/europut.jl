struct EuropeanPut <: AbstractContract
    K::Float64
    T::Float64
end

function u(state::BlackScholesState, contract::EuropeanPut)
    @unpack s, r, σ = state
    @unpack K, T = contract

    d1 = (log(s/K) + (r+0.5*σ^2)*T) / (√T*σ);
    d2 = d1 - σ*√T

    u = K*Φ(-d2) - exp(r*T)*s*Φ(-d1)
end

function h(state::BlackScholesState, contract::EuropeanPut)
    @unpack s = state
    @unpack K = contract

    h = max(K-s, 0)
end
