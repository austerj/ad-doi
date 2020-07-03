function estimator(nsteps, npaths, model::Heston, contract::AbstractContract, rng::AbstractRNG)
    @unpack s₀, ν₀, r, κ, θ, ξ, ρ = model
    @unpack T = contract

    Δ = T/nsteps
    u₀ = u(0., s₀, ν₀, model, contract)
    A₀ = diffop(0., s₀, ν₀, model, contract)/2

    payoffs = Vector{Float64}(undef, npaths)
    doi = Vector{Float64}(undef, npaths)

    for j = 1:npaths
        t, s, ν = 0., s₀, ν₀
        A = A₀
        for i = 1:nsteps
            t += Δ
            s, ν = step(s, ν, Δ, model, rng)
            A += diffop(t, s, ν, model, contract)
        end
        # trapezoidal rule; diffop is zero at time T hence division not needed
        A *= Δ

        payoffs[j] = h(s, contract)
        doi[j] = A
    end

    discount = exp(-r*T)
    mc_estimate = discount * mean(payoffs)
    doi_estimate = discount * (u₀ + mean(doi))

    mc_estimate, doi_estimate
end
