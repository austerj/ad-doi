using RandomNumbers.Xorshifts
using Parameters
using MuladdMacro
using ForwardDiff: Dual

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractContract end

using Distributions: cdf, pdf, Normal
Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

struct BlackScholes <: AbstractModel end
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("blackscholes/strangle.jl")
include("blackscholes/powercall.jl")
include("blackscholes/powerbinaryput.jl")

include("heston.jl")
include("diffop.jl")
include("estimator.jl")
