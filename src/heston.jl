@with_kw struct Heston <: AbstractModel
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

@muladd function step(s, ν, Δ, model::Heston, rng::AbstractRNG)
    @unpack r, κ, θ, ξ, ρ = model

    ΔW = √Δ*randn(rng)
    ΔZ = ρ*ΔW + √(1-ρ^2)*√Δ*randn(rng)

    # Euler-Maruyama predictors
    s̃ = s*(1 + Δ*r + √ν*ΔW)
    ν̃ = ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Drift-implicit correctors
    s += Δ/2*r*(s̃+s) + s*√ν*ΔW
    ν = max(ν + Δ*κ*(θ-(ν̃+ν)/2) + ξ*√ν*ΔZ, 1e-8)

    s, ν
end

function path(T, nsteps, npaths, model::Heston, rng::AbstractRNG)
    @unpack s₀, ν₀ = model

    Δ = T/nsteps

    s = Array{Float64,2}(undef, nsteps+1, npaths)
    ν = Array{Float64,2}(undef, nsteps+1, npaths)

    s[1,:] .= s₀
    ν[1,:] .= ν₀

    for j = 1:npaths
        s_step, ν_step = s₀, ν₀
        for i = 1:nsteps
            s_step, ν_step = step(s_step, ν_step, Δ, model, rng)
            s[i+1,j] = s_step
            ν[i+1,j] = ν_step
        end
    end

    s, ν
end

@muladd function u(t, s, ν, model::Heston, contract::AbstractContract, state::AbstractState)
    @unpack r, κ, θ, ξ, ρ = model
    @unpack T = contract
    τ = T-t

    if τ > 0
        # analytical solution to √(1/τ*∫ν̄ₜdt) from 0 to τ
        σ̄ = √((θ + (ν-θ)/κ*(1-exp(-κ*τ))/τ))
        u(t, s, σ̄, r, BlackScholes(), contract, state)
    else
        h(s, contract, state)
    end
end

# auto default state
function u(t, s, ν, model::Heston, contract::AbstractContract)
    u(t, s, ν, model::Heston, contract::AbstractContract, DefaultState())
end
