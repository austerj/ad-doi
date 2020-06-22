using Test

include("../base.jl")

# contract properties
T = 1
K = 90
n = 2

# heston properties
s₀ = 100
ν₀ = 0.16
r = 0.04
κ = 0.6
θ = 0.04
ξ = 0.2
ρ = -0.15

# black-scholes properties
σ = 0.2
r = 0.04

# terminal state
sT = 95

# tolerance for numerical comparisons
atol = 1e-6

@testset "Tests" begin
    @testset "Black-Scholes" begin include("blackscholes_test.jl") end
    @testset "Generalized Black-Scholes" begin include("generalizedblackscholes_test.jl") end
    @testset "Payoffs" begin include("payoff_test.jl") end
    @testset "Greeks" begin include("greeks_test.jl") end
end
