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

struct BlackScholes <: AbstractModel end
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("blackscholes/powercall.jl")
include("blackscholes/powerbinaryput.jl")

include("heston.jl")
include("diffop.jl")
include("estimator.jl")
