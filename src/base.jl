using RandomNumbers.Xorshifts
using Parameters
using MuladdMacro
using ForwardDiff: Dual
using BenchmarkTools
using SpecialFunctions: zeta

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractContract end
abstract type AbstractState end

using Distributions: cdf, pdf, Normal, std
Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

# bias correction factor
const β = -zeta(.5) / √(2*π)

include("blackscholes/model.jl")
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("blackscholes/powercall.jl")
include("blackscholes/powerbinaryput.jl")
include("blackscholes/lookback.jl")
include("blackscholes/downout_barrier_call.jl")
include("blackscholes/discrete_lookback.jl")

include("heston.jl")
include("diffop.jl")
include("estimator.jl")
