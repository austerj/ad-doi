include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes
using LaTeXStrings

pyplot()

rng = Xoroshiro128Plus(1)
fine_multiplier = 200
nsteps = 40
nsteps_fine = nsteps*fine_multiplier

T = 1.
t = range(0, stop=T, length=nsteps+1)
t_fine = range(0, stop=T, length=nsteps_fine+1)
Δ = T/nsteps
Δ_fine = T/nsteps_fine

theme(
    :default,
    size=(800, 300),
    xminorticks=4,
    yminorticks=3,
    xgrid=false,
    ygrid=false,
    legend=false,
    linewidth=1.5
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]
alpha = 0.8

# numerical schemes
function euler_maruyama_step(s, ν, Δ, model::Heston, ΔW, ΔZ)
    @unpack r, κ, θ, ξ = model
    s = s*(1 + Δ*r + √ν*ΔW)
    ν = max(ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ, 1e-8)
    s, ν
end

function predictor_corrector_step(s, ν, Δ, model::Heston, ΔW, ΔZ)
    @unpack r, κ, θ, ξ = model

    # Euler-Maruyama predictors
    s̃ = s*(1 + Δ*r + √ν*ΔW)
    ν̃ = ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Implicit ν corrector
    ν_drift = κ*(2*θ-(ν̃+ν)) - 0.5*ξ^2
    ν_vol = ξ*(√ν̃ + √ν)
    ν̂ = max(ν + .5*(ν_drift*Δ + ν_vol*ΔZ), 1e-8)

    # Implicit s corrector
    s_drift = s̃*(r-0.5*ν̃) + s*(r-0.5*ν)
    s_vol = s̃*√ν̃ + s*√ν
    ŝ = s + .5*(s_drift*Δ + s_vol*ΔW)

    ŝ, ν̂
end

# initialize paths
em = Array{Float64}(undef, nsteps+1, 2)
em_fine = Array{Float64}(undef, nsteps_fine+1, 2)
pc = Array{Float64}(undef, nsteps+1, 2)

em[1,1] = s₀
em[1,2] = ν₀
em_fine[1,1] = s₀
em_fine[1,2] = ν₀
pc[1,1] = s₀
pc[1,2] = ν₀

# plot
for i = 1:nsteps
    ΔW = 0
    ΔZ = 0
    offset = (i-1)*fine_multiplier
    for j in 1:fine_multiplier
        ΔW_fine = √Δ_fine*randn(rng)
        ΔZ_fine = ρ*ΔW_fine + √(1-ρ^2)*√Δ_fine*randn(rng)
        em_fine[offset+j+1,:] .= euler_maruyama_step(em_fine[offset+j,1], em_fine[offset+j,2], Δ_fine, heston, ΔW_fine, ΔZ_fine)
        ΔW += ΔW_fine
        ΔZ += ΔZ_fine
    end
    em[i+1,:] .= euler_maruyama_step(em[i,1], em[i,2], Δ, heston, ΔW, ΔZ)
    pc[i+1,:] .= predictor_corrector_step(pc[i,1], pc[i,2], Δ, heston, ΔW, ΔZ)
end

plot(t_fine, em_fine[:,1], linestyle=:dash, color="grey", label="euler-maruyama (fine)", alpha=0.5, linewidth=1)
plot!(t, em[:,1], color=colors[1], label="euler-maruyama", alpha=alpha, legend=true)
plot!(t, pc[:,1], color=colors[2], label="predictor corrector", alpha=alpha)
plot!(xlabel=L"t")

savefig("PredictorCorrector.pdf")
