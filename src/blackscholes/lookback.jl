struct FloatingLookbackPut <: AbstractContract
    T::Float64
end

# track running maximum
mutable struct RunningMax <: AbstractState
    max::Float64
end

function state(model::AbstractModel, contract::FloatingLookbackPut)
    RunningMax(model.s₀)
end

function state!(t, s, ν, contract::AbstractContract, state::RunningMax)
    if s > state.max
        state.max = s
    end
end

@muladd function u(t, s, σ, r, model::BlackScholes, contract::FloatingLookbackPut, state::RunningMax)
    @unpack T = contract
    @unpack max = state
    τ = T-t

    d1 = (log(s/max) + (r+0.5*σ^2)*τ) / (√τ*σ)
    d2 = d1 - σ*√τ
    d3 = d1 - 2*r*√τ/σ

    -exp(r*τ)*s*Φ(-d1) + max*Φ(-d2) + 0.5*s*σ^2/r*(exp(r*τ)*Φ(d1)-(max/s)^(2*r/σ^2)*Φ(d3))
end

function h(s, contract::FloatingLookbackPut, state::RunningMax)
    state.max - s
end
