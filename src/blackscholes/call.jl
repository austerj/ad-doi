struct EuropeanCall <: AbstractContract
    T::Float64
    K::Float64
end

function u(state::BlackScholesState, contract::EuropeanCall)
    @unpack s, σ, r = state
    @unpack T, K = contract

    d1 = (log(s/K) + (r+0.5*σ^2)*T) / (√T*σ);
    d2 = d1 - σ*√T

    u = exp(r*T)*s*Φ(d1) - K*Φ(d2)
end

function h(state::BlackScholesState, contract::EuropeanCall)
    @unpack s = state
    @unpack K = contract

    h = max(s-K, 0)
end
