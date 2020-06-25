# comparison to ensure correct implementation of value functions
const BLACKSCHOLES_PUT = 6.24902542413863
const BLACKSCHOLES_CALL = 10.3301028433775
const BLACKSCHOLES_POW_CALL = undef
const BLACKSCHOLES_BIN_CASH_CALL = undef

@test u(0., s₀, σ, r, BlackScholes(), EuropeanPut(T,K)) ≈ BLACKSCHOLES_PUT atol=atol
@test u(0., s₀, σ, r, BlackScholes(), EuropeanCall(T,K)) ≈ BLACKSCHOLES_CALL atol=atol
# @test u(0., s₀, σ, r, BlackScholes(), PowerCall(T,K,n)) ≈ BLACKSCHOLES_POWCALL atol=atol
# @test u(0., s₀, σ, r, BlackScholes(), BinaryCashCall(T,K)) ≈ BLACKSCHOLES_BIN atol=atol
