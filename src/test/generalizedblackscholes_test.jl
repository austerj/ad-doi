# ensure correct implementation of heston propagation to black-scholes
const GENERALIZEDBLACKSCHOLES_PUT =  12.6597977999531
const GENERALIZEDBLACKSCHOLES_CALL = 16.7408752191919
const GENERALIZEDBLACKSCHOLES_POW_CALL = undef
const GENERALIZEDBLACKSCHOLES_BIN_CASH_CALL = undef

@test u(0, s₀, ν₀, heston, EuropeanPut(T,K)) ≈ GENERALIZEDBLACKSCHOLES_PUT atol=atol
@test u(0, s₀, ν₀, heston, EuropeanCall(T,K)) ≈ GENERALIZEDBLACKSCHOLES_CALL atol=atol
# @test u(0, s₀, ν₀, heston, PowerCall(T,K,n)) ≈ GENERALIZEDBLACKSCHOLES_POW atol=atol
# @test u(0, s₀, ν₀, heston, Binary(T,K)) ≈ GENERALIZEDBLACKSCHOLES_BIN atol=atol
