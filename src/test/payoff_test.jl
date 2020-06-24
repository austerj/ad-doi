# payoffs to test that terminal prices = payoff
const PAYOFF_PUT = max(K-sT, 0)
const PAYOFF_CALL = max(sT-K, 0)
const PAYOFF_POW_CALL = max(sT^n-K, 0)
const PAYOFF_BIN_CASH_CALL = sT>K
const PAYOFF_POW_BIN_PUT = (sT<K) * sT^n

@test h(sT, EuropeanPut(T,K)) ≈ PAYOFF_PUT atol=atol
@test h(sT, EuropeanCall(T,K)) ≈ PAYOFF_CALL atol=atol
# @test h(sT, PowerCall(T,K,n)) ≈ PAYOFF_POWCALL atol=atol
# @test h(sT, BinaryCashCall(T,K)) PAYOFF_BIN atol=atol
# @test h(sT, PowerBinaryPut(T,K,n)) PAYOFF_POW_BIN_PUT atol=atol

# generalized black-scholes
@test u(T, sT, ν₀, heston, EuropeanPut(T,K)) ≈ PAYOFF_PUT atol=atol
@test u(T, sT, ν₀, heston, EuropeanCall(T,K)) ≈ PAYOFF_CALL atol=atol
# @test u(T, sT, ν₀, heston, PowerCall(T,K,n)) ≈ PAYOFF_POW atol=atol
# @test u(T, sT, ν₀, heston, Binary(T,K)) ≈ PAYOFF_BIN atol=atol
# @test u(T, sT, ν₀, heston, PowerBinaryCall(T,K)) ≈ PAYOFF_POW_BIN_PUT atol=atol
