# comparison to ensure correct implementation of value functions
const BLACKSCHOLES_PUT = 6.24902542413863
const BLACKSCHOLES_CALL = 10.3301028433775
const BLACKSCHOLES_POW_CALL = 11174.9685157938
const BLACKSCHOLES_POW_BIN_PUT = 3478.751035076

@test u(0., s₀, σ, r, BlackScholes(), EuropeanPut(T,K)) ≈ BLACKSCHOLES_PUT atol=atol
@test u(0., s₀, σ, r, BlackScholes(), EuropeanCall(T,K)) ≈ BLACKSCHOLES_CALL atol=atol
@test u(0., s₀, σ, r, BlackScholes(), PowerCall(T,K,n)) ≈ BLACKSCHOLES_POW_CALL atol=atol
@test u(0., s₀, σ, r, BlackScholes(), PowerBinaryPut(T,K,n)) ≈ BLACKSCHOLES_POW_BIN_PUT atol=atol
