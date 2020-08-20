@muladd function diffop(t, s, ν, model::Heston, contract::AbstractContract, state::AbstractState)::Float64
    @unpack r, κ, θ, ξ, ρ = model

    ∂sν, ∂νν = sensitivities(t, s, ν, model, contract, state)
    ξ*ν*(ρ*s*∂sν + 0.5*ξ*∂νν)
end

function sensitivities(t, s, ν, model, contract, state)::Tuple{Float64,Float64}
    zs = Dual(Dual(s, 1., 0.), Dual(0., 0., 0.))
    zν = Dual(Dual(ν, 0., 1.), Dual(1., 0., 0.))
    uz = u(t, zs, zν, model, contract, state)

    ∂sν = uz.partials[1].partials[1]  # vanna
    ∂νν = uz.partials[1].partials[2]  # vomma

    ∂sν, ∂νν
end
