@muladd function diffop(t, s, ν, model::MultivariateHeston, contract::AbstractContract)::Float64
    @unpack ξ, ρ = model

    Hsν, Hνν = sensitivities(t, s, ν, model, contract)

    Ds = Diagonal(s)
    Dξ = Diagonal(ξ)

    d = length(s)
    ρsν = ρ[d+1:end, 1:d]
    ρνν = ρ[d+1:end, d+1:end]

    transpose(.√ν) * (Ds*(ρsν.*Hsν) + 0.5*Dξ*(ρνν.*Hνν)) * Dξ*.√ν
end

function sensitivities(t, s, ν, model::MultivariateHeston, contract::AbstractContract)
    X = [s; ν]
    H = hessian(x -> u(t, x, model, contract), X)

    d = length(s)
    Hsν = H[d+1:end, 1:d]  # vanna
    Hνν = H[d+1:end, d+1:end]  # vomma

    Hsν, Hνν
end
