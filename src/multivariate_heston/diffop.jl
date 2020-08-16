@muladd function diffop(t, s, ν, model::MultivariateHeston, contract::AbstractContract)::Float64
    @unpack ξ, ρ, cholesky, ndims = model

    Hsν, Hνν = sensitivities(t, s, ν, model, contract)

    # matrix implementation - significantly slower than direct sum
    # Ds = Diagonal(s)
    # Dξ = Diagonal(ξ)

    # d = length(s)
    # ρsν = @view ρ[d+1:end, 1:d]
    # ρνν = @view ρ[d+1:end, d+1:end]

    # transpose(.√ν) * (Ds*(ρsν.*Hsν) + 0.5*Dξ*(ρνν.*Hνν)) * Dξ*.√ν

    ddiffop = 0
    for i=1:ndims, j=1:ndims
        factor = ξ[i]*√(ν[i]*ν[j])
        term_sν = s[i]*Hsν[i,j]
        term_νν = .5*ξ[i]*Hνν[i,j]

        sum_k = 0
        for k=1:ndims
             sum_k += term_sν*cholesky[k+i,j] + term_νν*cholesky[k+i,k+j]
        end

        ddiffop += factor * sum_k
    end

    ddiffop
end

function sensitivities(t, s, ν, model::MultivariateHeston, contract::AbstractContract)
    @unpack ndims = model

    H = hessian(x -> u(t, x, model, contract), [s; ν])

    Hsν = H[ndims+1:end, 1:ndims]  # vanna
    Hνν = H[ndims+1:end, ndims+1:end]  # vomma

    Hsν, Hνν
end
