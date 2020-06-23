# comparison to ensure correct implementation of value functions
const BLACKSCHOLES_PUT = 6.24902542413863
const BLACKSCHOLES_CALL = 10.3301028433775
const BLACKSCHOLES_POW_CALL = undef
const BLACKSCHOLES_BIN_CASH_CALL = undef

state = BlackScholesState(0, s₀, σ, r)

@test u(state, EuropeanPut(T,K)) ≈ BLACKSCHOLES_PUT atol=atol
@test u(state, EuropeanCall(T,K)) ≈ BLACKSCHOLES_CALL atol=atol
# @test u(state, PowerCall(T,K,n)) ≈ BLACKSCHOLES_POWCALL atol=atol
# @test u(state, BinaryCashCall(T,K)) ≈ BLACKSCHOLES_BIN atol=atol
