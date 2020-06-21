using Test

include("../base.jl")

# contract properties
const T = 1
const K = 90
const n = 2

# heston properties
const s₀ = 100
const ν₀ = 0.16
const r = 0.04
const κ = 0.6
const θ = 0.04
const ξ = 0.2
const ρ = -0.15

# black-scholes properties
const σ = 0.2
const r = 0.04

# terminal state
const sT = 95

# tolerance for numerical comparisons
const atol = 1e-6

@testset "Tests" begin
    @testset "Black-Scholes" begin include("blackscholes_test.jl") end
    @testset "Generalized Black-Scholes" begin include("generalizedblackscholes_test.jl") end
    @testset "Payoffs" begin include("payoff_test.jl") end
    @testset "Greeks" begin include("greeks_test.jl") end
end
