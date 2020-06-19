using RandomNumbers.Xorshifts
using BenchmarkTools
using Parameters
using MuladdMacro

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractPath end
abstract type AbstractContract end
abstract type AbstractState end

using Distributions: cdf, Normal
Φ(x) = cdf(Normal(), x)

include("blackscholes/blackscholes.jl")
include("blackscholes/europut.jl")
include("heston.jl")
