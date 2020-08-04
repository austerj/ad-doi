using RandomNumbers.Xorshifts
using Parameters
using MuladdMacro
using LinearAlgebra
using ForwardDiff: Dual, hessian

using RandomNumbers: AbstractRNG
abstract type AbstractModel end
abstract type AbstractContract end

using Distributions: cdf, pdf, Normal, std
Φ(x) = cdf(Normal(), x)
φ(x) = pdf(Normal(), x)

struct BlackScholes <: AbstractModel end
include("blackscholes/put.jl")
include("blackscholes/call.jl")
include("blackscholes/strangle.jl")
include("blackscholes/powercall.jl")
include("blackscholes/powerbinaryput.jl")
include("blackscholes/smoothedpowerbinaryput.jl")

include("heston/heston.jl")
include("heston/diffop.jl")
include("heston/estimator.jl")

include("multivariate_heston/multivariate_heston.jl")
include("multivariate_heston/diffop.jl")
include("multivariate_heston/estimator.jl")
