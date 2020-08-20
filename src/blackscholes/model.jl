struct BlackScholes <: AbstractModel end
struct DefaultState <: AbstractState end

# default state for contracts without path dependence
function state(model::AbstractModel, contract::AbstractContract)
    DefaultState()
end

function state!(t, s, ν, contract::AbstractContract, state::DefaultState) end

# call value function for contracts using DefaultState
function u(t, s, σ, r, model::BlackScholes, contract::AbstractContract, state::DefaultState)
    u(t, s, σ, r, model, contract)
end

# call payoff functions for contracts using DefaultState
function h(s, contract::AbstractContract, state::DefaultState)
    h(s, contract)
end
