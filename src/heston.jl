using Parameters
using MuladdMacro
using RandomNumbers: AbstractRNG

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

mutable struct HestonState
    s::Float64
    ν::Float64
end

@muladd function step!(state::HestonState, Δ, model::Heston, rng::AbstractRNG)
    @unpack s, ν = state
    @unpack r, κ, θ, ξ, ρ = model

    ΔW = √Δ*randn(rng)
    ΔZ = √(1-ρ^2)*√Δ*randn(rng)

    # Euler-Maruyama predictors
    s̃ = s*(1 + Δ*r + √ν*ΔW)
    ν̃ = ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Drift-implicit correctors
    s += Δ/2*r*(s̃+s) + s*√ν*ΔW
    ν = max(ν + Δ*κ*(θ-(ν̃+ν)/2) + ξ*√ν*ΔZ, 1e-14)

    state.s, state.ν = s, ν
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

    for j=1:npaths
        state = HestonState(s₀, ν₀)
        for i=1:nsteps
            step!(state, Δ, model, rng)
            s[i+1,j] = state.s
            ν[i+1,j] = state.ν
        end
    end

    s, ν
end
