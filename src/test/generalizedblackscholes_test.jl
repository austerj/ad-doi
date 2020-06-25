# ensure correct implementation of heston propagation to black-scholes
const GENERALIZEDBLACKSCHOLES_PUT =  12.6597977999531
const GENERALIZEDBLACKSCHOLES_CALL = 16.7408752191919
const GENERALIZEDBLACKSCHOLES_POW_CALL = 12239.7130649947
const GENERALIZEDBLACKSCHOLES_POW_BIN_PUT = 3173.12226966096

@test u(0., s₀, ν₀, heston, EuropeanPut(T,K)) ≈ GENERALIZEDBLACKSCHOLES_PUT atol=atol
@test u(0., s₀, ν₀, heston, EuropeanCall(T,K)) ≈ GENERALIZEDBLACKSCHOLES_CALL atol=atol
@test u(0., s₀, ν₀, heston, PowerCall(T,K,n)) ≈ GENERALIZEDBLACKSCHOLES_POW_CALL atol=atol
@test u(0., s₀, ν₀, heston, PowerBinaryPut(T,K,n)) ≈ GENERALIZEDBLACKSCHOLES_POW_BIN_PUT atol=atol
