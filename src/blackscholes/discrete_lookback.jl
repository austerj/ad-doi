struct DiscreteFloatingLookbackPut <: AbstractContract
    T::Float64
    nsteps::Int32
end

# track running maximum
mutable struct RunningMax <: AbstractState
    max::Float64
end

function state(model::AbstractModel, contract::DiscreteFloatingLookbackPut)
    RunningMax(model.s₀)
end

function state!(t, s, ν, contract::DiscreteFloatingLookbackPut, state::RunningMax)
    if s > state.max
        state.max = s
    end
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::DiscreteFloatingLookbackPut, state::RunningMax)
    @unpack T, nsteps = contract
    @unpack max = state
    τ = T-t

    # bias correction for discretely-monitored price
    bias_exponent = β*σ*√(T/nsteps)
    
    # bias correct running max
    mod_max = max*exp(bias_exponent)

    d1 = (log(s/mod_max) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ
    d3 = d1 - 2*r*√τ/σ

    value = -exp(r*τ)*s*Φ(-d1) + mod_max*Φ(-d2) + 0.5*s*σ^2/r*(exp(r*τ)*Φ(d1)-(mod_max/s)^(2*r/σ^2)*Φ(d3))

    # bias corrected value function
    exp(-bias_exponent)*value + s*exp(r*τ)*(exp(-bias_exponent)-1)
end

function h(s, contract::DiscreteFloatingLookbackPut, state::RunningMax)
    state.max - s
end
