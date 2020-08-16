@with_kw struct MultivariateHeston <: AbstractModel
    s₀::Vector{Float64}
    ν₀::Vector{Float64}
    r::Float64
    κ::Vector{Float64}
    θ::Vector{Float64}
    ξ::Vector{Float64}

    # ρ is the correlation matrix of the form
    #    W1 W2 . . . Zd
    # W1 1 
    # W2    1
    # .        .
    # .          .
    # .            .
    # Zd             1
    ρ::Array{Float64,2}
    cholesky::LowerTriangular

    ndims::Integer

    function MultivariateHeston(s₀, ν₀, r, κ, θ, ξ, ρ)
        ndims = length(s₀)
        if any([length(x) != ndims for x in [s₀, ν₀, κ, θ, ξ]]) || size(ρ) != (2*ndims, 2*ndims)
            error("Dimension mismatch")
        elseif any(2*κ.*θ <= ξ.^2)
            error("Feller condition not satisfied")
        else
            new(s₀, ν₀, r, κ, θ, ξ, Symmetric(ρ, :U), cholesky(ρ).L, ndims)
        end
    end
end

@muladd function step(s, ν, Δ, model::MultivariateHeston, rng::AbstractRNG)::Tuple{Vector{Float64},Vector{Float64}}
    @unpack r, κ, θ, ξ, cholesky, ndims = model

    ΔB = cholesky * √Δ*randn(rng, 2*ndims)
    ΔW = @view ΔB[1:ndims]
    ΔZ = @view ΔB[ndims+1:end]

    # Euler-Maruyama predictors
    s̃ = @. s*(1 + Δ*r + √ν*ΔW)
    ν̃ = @. ν + Δ*κ*(θ-ν) + ξ*√ν*ΔZ

    # Drift-implicit correctors
    s += @. Δ/2*r*(s̃+s) + s*√ν*ΔW
    ν = @. max(ν + Δ*κ*(θ-(ν̃+ν)/2) + ξ*√ν*ΔZ, 1e-8)

    s, ν
end

@muladd function u(t, x, model::MultivariateHeston, contract::AbstractContract)
    @unpack r, κ, θ, ξ, ρ, ndims = model
    @unpack T = contract
    τ = T-t
    
    s = x[1:ndims]

    if τ > 0
        ν = x[ndims+1:2*ndims]

        # analytical solution to √(1/τ*∫ν̄ₜdt) from 0 to τ
        σ̄ = @. √((θ + (ν-θ)/κ*(1-exp(-κ*τ))/τ))
        u(t, s, σ̄, r, ρ[1:ndims, 1:ndims], BlackScholes(), contract)
    else
        h(s, contract)
    end
end
