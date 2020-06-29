@muladd function diffop(t, s, ν, model::Heston, contract::AbstractContract)::Float64
    @unpack r, κ, θ, ξ, ρ = model

    ∂sν, ∂νν = greeks(t, s, ν, model, contract)
    ξ*ν*(ρ*s*∂sν + 0.5*ξ*∂νν)
end

function greeks(t, s, ν, model::Heston, contract::AbstractContract)::Tuple{Float64,Float64}
    Ds = Dual(Dual(s, 1., 0.), Dual(0., 0., 0.))
    Dν = Dual(Dual(ν, 0., 1.), Dual(1., 0., 0.))
    greeks = u(t, Ds, Dν, model, contract)

    ∂sν = greeks.partials[1].partials[1]
    ∂νν = greeks.partials[1].partials[2]

    ∂sν, ∂νν
end
