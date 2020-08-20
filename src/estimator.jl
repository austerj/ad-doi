function estimator(nsteps, npaths, model::Heston, contract::AbstractContract, rng::AbstractRNG)
    @unpack s₀, ν₀, r, κ, θ, ξ, ρ = model
    @unpack T = contract

    Δ = T/nsteps
    t = Δ:Δ:T
    initial_state = state(model, contract)
    u₀ = u(0., s₀, ν₀, model, contract, initial_state)
    A₀ = diffop(0., s₀, ν₀, model, contract, initial_state)/2

    payoffs = Vector{Float64}(undef, npaths)
    doi = Vector{Float64}(undef, npaths)

    for j = 1:npaths
        s, ν = s₀, ν₀
        A = A₀
        path_state = deepcopy(initial_state)
        for i = 1:nsteps
            s, ν = step(s, ν, Δ, model, rng)
            state!(t, s, ν, contract, path_state)
            A += diffop(t[i], s, ν, model, contract, path_state)
        end
        # trapezoidal rule; diffop is zero at time T hence division not needed
        A *= Δ

        payoffs[j] = h(s, contract, path_state)
        doi[j] = A
    end

    discount = exp(-r*T)
    mc_estimate = discount * mean(payoffs)
    doi_estimate = discount * (u₀ + mean(doi))

    mc_estimate, doi_estimate
end

# set threads with "export JULIA_NUM_THREADS=$(nproc)" outside REPL
function parallel_estimator(nsteps, npaths, model::Heston, contract::AbstractContract)
    @unpack s₀, ν₀, r = model
    @unpack T = contract

    Δ = T/nsteps
    t = Δ:Δ:T
    initial_state = state(model, contract)
    u₀ = u(0., s₀, ν₀, model, contract, initial_state)
    A₀ = Δ*diffop(0., s₀, ν₀, model, contract, initial_state)/2

    nthreads = Threads.nthreads()
    npath_atomic = Threads.Atomic{Int}(0)
    payoff_atomic = Threads.Atomic{Float64}(0)
    A_atomic = Threads.Atomic{Float64}(0)

    # estimates pass two-sample Kolmogorov-Smirnov test despite no randjump
    rngs = [Xoroshiro128Plus(rand(UInt128)) for i = 1:nthreads]

    Threads.@threads for i in 1:Threads.nthreads()
        while true
            npath_thread = Threads.atomic_add!(npath_atomic, 1)

            if npath_thread < npaths
                payoff, A = sample(t, Δ, nsteps, model, contract, initial_state, rngs[Threads.threadid()])
                Threads.atomic_add!(A_atomic, A) 
                Threads.atomic_add!(payoff_atomic, payoff) 
            else
                break
            end
        end
    end

    discount = exp(-r*T)
    mc_estimate = discount * payoff_atomic[]/npaths
    doi_estimate = discount * (u₀ + A₀ + A_atomic[]/npaths)

    mc_estimate, doi_estimate
end

function sample(t, Δ, nsteps, model::Heston, contract::AbstractContract, initial_state::AbstractState, rng::AbstractRNG)
    @unpack s₀, ν₀, r, κ, θ, ξ, ρ = model

    s, ν = s₀, ν₀
    A = 0
    path_state = deepcopy(initial_state)
    @inbounds for i = 1:nsteps
        s, ν = step(s, ν, Δ, model, rng)
        state!(t, s, ν, contract, path_state)
        A += diffop(t[i], s, ν, model, contract, path_state)
    end
    # trapezoidal rule; diffop is zero at time T hence division not needed
    A *= Δ
    payoff = h(s, contract, path_state)

    payoff, A
end

# mc-only estimator
function mc_parallel_estimator(nsteps, npaths, model::Heston, contract::AbstractContract)
    @unpack s₀, ν₀, r = model
    @unpack T = contract

    Δ = T/nsteps
    t = Δ:Δ:T
    initial_state = state(model, contract)

    nthreads = Threads.nthreads()
    npath_atomic = Threads.Atomic{Int}(0)
    payoff_atomic = Threads.Atomic{Float64}(0)

    # estimates pass two-sample Kolmogorov-Smirnov test despite no randjump
    rngs = [Xoroshiro128Plus(rand(UInt128)) for i = 1:nthreads]

    Threads.@threads for i in 1:Threads.nthreads()
        while true
            npath_thread = Threads.atomic_add!(npath_atomic, 1)

            if npath_thread < npaths
                payoff = mc_sample(t, Δ, nsteps, model, contract, initial_state, rngs[Threads.threadid()])
                Threads.atomic_add!(payoff_atomic, payoff) 
            else
                break
            end
        end
    end

    discount = exp(-r*T)
    mc_estimate = discount * payoff_atomic[]/npaths

    mc_estimate
end

function mc_sample(t, Δ, nsteps, model::Heston, contract::AbstractContract, initial_state::AbstractState, rng::AbstractRNG)
    @unpack s₀, ν₀, r, κ, θ, ξ, ρ = model

    s, ν = s₀, ν₀
    path_state = deepcopy(initial_state)
    @inbounds for i = 1:nsteps
        s, ν = step(s, ν, Δ, model, rng)
        state!(t, s, ν, contract, path_state)
    end

    h(s, contract, path_state)
end
