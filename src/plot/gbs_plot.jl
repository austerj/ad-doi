include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes
using LaTeXStrings

pyplot()

rng = Xoroshiro128Plus(1)
nsteps = 400
npaths = 1

T = 5.
t = range(0, stop=T, length=nsteps+1)
s, ν = path(T, nsteps, npaths, heston, rng)

theme(
    :default,
    size=(800, 250),
    xminorticks=4,
    yminorticks=3,
    xgrid=false,
    ygrid=false,
    legend=false
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]
stops = [Int(ceil(nsteps/k)) for k in [7, 3, 2]]
lwidth = 1
alpha = 0.8

p = plot()
for j = 1:npaths
    stop = stops[j]
    σ̄ = @. (θ + (ν[stop, j]-θ)*exp(-κ*(t[stop:end]-t[stop])))

    plot!(t[1:stop], ν[1:stop, j], color=colors[1], label=L"\nu", alpha=alpha)
    plot!(t[stop:end], ν[stop:end, j], color=colors[2], label="", alpha=0.5)
    plot!(t[stop:end], σ̄, linestyle=:dash, color=colors[1], label=L"\bar\nu", alpha=alpha)
end
plot!(xticks=false)

savefig("VarianceProcessApproximation.pdf")
