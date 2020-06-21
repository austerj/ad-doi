using Test

include("../base.jl")

# contract properties
const global T = 1
const global K = 90
const global n = 2

# tolerance for numerical comparisons
const global atol = 1e-6

@testset begin include("blackscholes_test.jl") end
@testset begin include("generalizedblackscholes_test.jl") end
