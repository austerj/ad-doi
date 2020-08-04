@muladd function diffop(t, s, ν, model::MultivariateHeston, contract::AbstractContract)::Float64
    @unpack ξ, ρ = model

    Hsν, Hνν = sensitivities(t, s, ν, model, contract)

    Ds = Diagonal(s)
    Dξ = Diagonal(ξ)

    d = length(s)
    ρsν = @view ρ[d+1:end, 1:d]
    ρνν = @view ρ[d+1:end, d+1:end]

    transpose(.√ν) * (Ds*(ρsν.*Hsν) + 0.5*Dξ*(ρνν.*Hνν)) * Dξ*.√ν
end

function sensitivities(t, s, ν, model::MultivariateHeston, contract::AbstractContract)
    @unpack ndims = model

    H = hessian(x -> u(t, x, model, contract), [s; ν])

    Hsν = H[ndims+1:end, 1:ndims]  # vanna
    Hνν = H[ndims+1:end, ndims+1:end]  # vomma

    Hsν, Hνν
end
