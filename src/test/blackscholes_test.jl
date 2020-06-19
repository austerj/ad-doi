# comparison to ensure correct implementation of value functions
const BLACKSCHOLES_PUT = 2.63478909638165
const BLACKSCHOLES_CALL = 16.7158665156205
const BLACKSCHOLES_POW_CALL = undef
const BLACKSCHOLES_BIN_CASH_CALL = undef
atol = 1e-6

# state
const s = 100
const σ = 0.2
const r = 0.04
state = BlackScholesState(s, σ, r)

# contract properties
const T = 1
const K = 90
const n = 2

@test u(state, EuropeanPut(T,K)) ≈ BLACKSCHOLES_PUT atol=atol
@test u(state, EuropeanCall(T,K)) ≈ BLACKSCHOLES_CALL atol=atol
# @test u(state, PowerCall(T,K,n)) ≈ BLACKSCHOLES_POWCALL atol=atol
# @test u(state, BinaryCashCall(T,K)) ≈ BLACKSCHOLES_BIN atol=atol
