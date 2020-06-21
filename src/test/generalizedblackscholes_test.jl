# comparison to ensure correct implementation of HestonPath type promotion
const GENERALIZEDBLACKSCHOLES_PUT =  7.95293256258722
const GENERALIZEDBLACKSCHOLES_CALL = 22.034009981826
const GENERALIZEDBLACKSCHOLES_POW_CALL = undef
const GENERALIZEDBLACKSCHOLES_BIN_CASH_CALL = undef

# heston
heston = Heston(s₀, ν₀, r, κ, θ, ξ, ρ)

@test u(0, s₀, ν₀, heston, EuropeanPut(T,K)) ≈ GENERALIZEDBLACKSCHOLES_PUT atol=atol
@test u(0, s₀, ν₀, heston, EuropeanCall(T,K)) ≈ GENERALIZEDBLACKSCHOLES_CALL atol=atol
# @test u(0, s₀, ν₀, heston, PowerCall(T,K,n)) ≈ GENERALIZEDBLACKSCHOLES_POW atol=atol
# @test u(0, s₀, ν₀, heston, Binary(T,K)) ≈ GENERALIZEDBLACKSCHOLES_BIN atol=atol
