# payoffs to test that terminal prices = payoff
const PAYOFF_PUT = max(K-sT, 0)
const PAYOFF_CALL = max(sT-K, 0)
const PAYOFF_POW_CALL = max(sT^n-K, 0)
const PAYOFF_POW_BIN_PUT = (sT<K) * sT^n

@test h(sT, EuropeanPut(T,K)) ≈ PAYOFF_PUT atol=atol
@test h(sT, EuropeanCall(T,K)) ≈ PAYOFF_CALL atol=atol
@test h(sT, PowerCall(T,K,n)) ≈ PAYOFF_POW_CALL atol=atol
@test h(sT, PowerBinaryPut(T,K,n)) ≈ PAYOFF_POW_BIN_PUT atol=atol

# generalized black-scholes
@test u(T, sT, ν₀, heston, EuropeanPut(T,K)) ≈ PAYOFF_PUT atol=atol
@test u(T, sT, ν₀, heston, EuropeanCall(T,K)) ≈ PAYOFF_CALL atol=atol
@test u(T, sT, ν₀, heston, PowerCall(T,K,n)) ≈ PAYOFF_POW_CALL atol=atol
@test u(T, sT, ν₀, heston, PowerBinaryPut(T,K,n)) ≈ PAYOFF_POW_BIN_PUT atol=atol
