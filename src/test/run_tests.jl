using Test

include("../base.jl")
include("../params.jl")

# contracts
T = 1.
K = 100.
n = 2.

# black-scholes
σ = 0.2
r = 0.04

# terminal state
sT = 110.

# estimator
rng = Xorshift128Plus(1)
nsteps = 50
npaths = 5000

# parallel test
npaths_run = 50
nruns = 1000
significance = 0.001

atol = 1e-6

@testset "Tests" begin
    @testset "Black-Scholes" begin include("blackscholes_test.jl") end
    @testset "Generalized Black-Scholes" begin include("generalizedblackscholes_test.jl") end
    @testset "Payoffs" begin include("payoff_test.jl") end
    @testset "Sensitivities" begin include("sensitivities_test.jl") end
    @testset "Estimator" begin include("estimator_test.jl") end
    @testset "Parallel Estimator" begin include("parallel_test.jl") end
    @testset "Multidimensional Heston" begin include("multi_test.jl") end
end
