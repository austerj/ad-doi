# comparison to ensure correct implementation of HestonPath type promotion
const GENERALIZEDBLACKSCHOLES_PUT =  7.95293256258722
const GENERALIZEDBLACKSCHOLES_CALL = 22.034009981826
const GENERALIZEDBLACKSCHOLES_POW_CALL = undef
const GENERALIZEDBLACKSCHOLES_BIN_CASH_CALL = undef

# heston
const s₀ = 100
const ν₀ = 0.16
const r = 0.04
const κ = 0.6
const θ = 0.04
const ξ = 0.2
const ρ = -0.15
heston = Heston(s, ν₀, r, κ, θ, ξ, ρ)

@test u(0, s₀, ν₀, heston, EuropeanPut(T,K)) ≈ GENERALIZEDBLACKSCHOLES_PUT atol=atol
@test u(0, s₀, ν₀, heston, EuropeanCall(T,K)) ≈ GENERALIZEDBLACKSCHOLES_CALL atol=atol
# @test u(path, PowerCall(r,T,K,n)) ≈ GENERALIZEDBLACKSCHOLES_POW atol=atol
# @test u(path, Binary(r,T)) ≈ GENERALIZEDBLACKSCHOLES_BIN atol=atol
