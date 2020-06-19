using Test

include("../base.jl")

# state
s = 100
σ = 0.2
r = 0.04
const global state = BlackScholesState(s, σ, r)

# contract properties
const global T = 1
const global K = 90
const global n = 2

# tolerance for numerical comparisons
const global atol = 1e-6

begin include("blackscholes_test.jl") end
