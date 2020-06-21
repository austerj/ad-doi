using RandomNumbers.Xorshifts
using BenchmarkTools
using Parameters
using MuladdMacro
using ForwardDiff: Dual

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractContract end
abstract type AbstractState end

using Distributions: cdf, Normal
Φ(x) = cdf(Normal(), x)

include("blackscholes/blackscholes.jl")
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("heston.jl")
