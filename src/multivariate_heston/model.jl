@with_kw struct MultivariateHeston <: AbstractModel
    s₀::Vector{Float64}
    ν₀::Vector{Float64}
    r::Float64
    κ::Vector{Float64}
    θ::Vector{Float64}
    ξ::Vector{Float64}

    # ρ is the correlation matrix of the form
    #    W1 W2 ... Z1 Z2
    # W1 1 
    # W2    1
    # .       .
    # .         .
    # .           .
    # Z1           1
    # Z2              1
    ρ::Array{Float64, 2}
    cholesky::Array{Float64, 2}

    ndims::Integer

    function MultivariateHeston(s₀, ν₀, r, κ, θ, ξ, ρ)
        any(2*κ.*θ <= ξ.^2) ? error("Feller condition not satisfied") : new(s₀, ν₀, r, κ, θ, ξ, ρ, cholesky(ρ).L, length(s₀))
    end
end

@muladd function step(s, ν, Δ, model::MultivariateHeston, rng::AbstractRNG)
    @unpack r, κ, θ, ξ, cholesky = model
    ndim = length(s)

    ΔB = cholesky * √Δ*randn(rng, 2*ndim)
    ΔW = @view ΔB[1:ndim]
    ΔZ = @view ΔB[ndim+1:end]

    # Euler-Maruyama predictors
    s̃ = @. s*(1 + Δ*r + √ν*ΔW)
    ν̃ = @. ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Drift-implicit correctors
    s += @. Δ/2*r*(s̃+s) + s*√ν*ΔW
    ν = @. max(ν + Δ*κ*(θ-(ν̃+ν)/2) + ξ*√ν*ΔZ, 1e-14)

    s, ν
end

@muladd function u(t, x, model::MultivariateHeston, contract::AbstractContract)
    @unpack r, κ, θ, ξ, ρ = model
    @unpack T = contract
    τ = T-t
    
    s = x[1]
    ν = x[2]

    if τ > 0
        # analytical solution to √(1/τ*∫ν̄ₜdt) from 0 to τ
        σ̄ = @. √((θ + (ν-θ)/κ*(1-exp(-κ*τ))/τ))
        u(t, s, σ̄, r, BlackScholes(), contract)
    else
        h(s, contract)
    end
end
