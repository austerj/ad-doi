# test that multidimensional heston with d=1 gives same estimate as normal heston
const NUMERICALINTEGRATION_PUT = 12.0190225503127

# define european put for multidimensional heston with d=1
struct _EuropeanPut <: AbstractContract
    T::Float64
    K::Float64
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::_EuropeanPut)
    @unpack T, K = contract
    τ = T-t

    s = s[1]
    σ = σ[1]

    d1 = (log(s/K) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ

    K*Φ(-d2) - exp(r*τ)*s*Φ(-d1)
end

function h(s, contract::_EuropeanPut)
    @unpack K = contract
    max(K-s[1], 0)
end

# parameter set H
_heston = MultivariateHeston([s₀], [ν₀], r, [κ], [θ], [ξ], [1 ρ; ρ 1])

# relative error tolerance of 0.1%
err_tol = 1e-3

# mean relative error ~1e-5 with nsteps=50, npaths=5000 based on 1000 samples for call
@test estimator(nsteps, npaths, _heston, _EuropeanPut(T,K), rng)[2] ≈ NUMERICALINTEGRATION_PUT atol=err_tol*NUMERICALINTEGRATION_PUT
