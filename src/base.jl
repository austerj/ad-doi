using RandomNumbers.Xorshifts
using Parameters
using MuladdMacro
using ForwardDiff: Dual

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractContract end
abstract type AbstractState end

using Distributions: cdf, pdf, Normal, std
Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

include("blackscholes/model.jl")
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("blackscholes/strangle.jl")
include("blackscholes/powercall.jl")
include("blackscholes/powerbinaryput.jl")
include("blackscholes/smoothedpowerbinaryput.jl")
include("blackscholes/lookback.jl")

include("heston.jl")
include("diffop.jl")
include("estimator.jl")
