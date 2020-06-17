using Parameters

struct Heston
    s₀::Float64
    ν₀::Float64
    r::Float64
    κ::Float64
    θ::Float64
    ξ::Float64
    ρ::Float64

    function Heston(s₀, ν₀, r, κ, θ, ξ, ρ)
        2*κ*θ <= ξ^2 ? error("Feller condition not satisfied") : new(s₀, ν₀, r, κ, θ, ξ, ρ)
    end
end

function step(s, ν, Δ, model::Heston, rng::AbstractRNG)::Tuple{Float64,Float64}
    @unpack r, κ, θ, ξ, ρ = model

    ΔW = √Δ*randn(rng)
    ΔZ = √(1-ρ^2)*√Δ*randn(rng)

    # Euler-Maruyama predictors
    s̃ = s*(1 + Δ*r + √ν*ΔW)
    ν̃ = ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Drift-implicit correctors
    s += Δ/2*r*(s̃+s) + s*√ν*ΔW
    ν = max(ν + Δ*κ*(θ-(ν̃+ν)/2) + ξ*√ν*ΔZ, 1e-14)

    s, ν
end

function path(T, nsteps, npaths, model::Heston, rng::AbstractRNG)
    @unpack s₀, ν₀ = model

    Δ = T/nsteps

    # preallocate arrays
    s = Array{Float64,2}(undef, nsteps+1, npaths)
    ν = Array{Float64,2}(undef, nsteps+1, npaths)

    # initialize first row
    s[1,:] .= s₀
    ν[1,:] .= ν₀

    for j=1:npaths, i=1:nsteps
        s_ij, ν_ij = step(s[i,j], ν[i,j], Δ, model, rng)
        s[i+1,j] = s_ij
        ν[i+1,j] = ν_ij
    end

    s, ν
end
