# test automatic differentiation against analytical solutions to undiscounted greeks
const BLACKSCHOLES_DELTA_PUT = -0.397681908481585
const BLACKSCHOLES_VEGA_PUT = 39.6952547477012
const BLACKSCHOLES_VANNA_PUT = -0.198476273738506
const BLACKSCHOLES_VOMMA_PUT = 5.95428821215518
const GENERALIZEDBLACKSCHOLES_DELTA_PUT = -0.40114740442307
const GENERALIZEDBLACKSCHOLES_VEGA_PUT = 41.463541599198
const GENERALIZEDBLACKSCHOLES_VANNA_PUT = 0.0799704029843022
const GENERALIZEDBLACKSCHOLES_VOMMA_PUT = -122.130268204048

# black-scholes greeks
Ds = Dual(Dual(s₀,1,0), Dual(0,0,0))
Dσ = Dual(Dual(σ,0,1), Dual(1,0,0))
bsgreeks = u(0., Ds, Dσ, r, BlackScholes(), EuropeanPut(T,K))

@test bsgreeks.value.partials[1] ≈ BLACKSCHOLES_DELTA_PUT atol=atol
@test bsgreeks.value.partials[2] ≈ BLACKSCHOLES_VEGA_PUT atol=atol
@test bsgreeks.partials[1].partials[1] ≈ BLACKSCHOLES_VANNA_PUT atol=atol
@test bsgreeks.partials[1].partials[2] ≈ BLACKSCHOLES_VOMMA_PUT atol=atol

# generalized black-scholes
Dν = Dual(Dual(ν₀,0,1), Dual(1,0,0))
gbsgreeks = u(0., Ds, Dν, heston, EuropeanPut(T,K))

@test gbsgreeks.value.partials[1] ≈ GENERALIZEDBLACKSCHOLES_DELTA_PUT atol=atol
@test gbsgreeks.value.partials[2] ≈ GENERALIZEDBLACKSCHOLES_VEGA_PUT atol=atol
@test gbsgreeks.partials[1].partials[1] ≈ GENERALIZEDBLACKSCHOLES_VANNA_PUT atol=atol
@test gbsgreeks.partials[1].partials[2] ≈ GENERALIZEDBLACKSCHOLES_VOMMA_PUT atol=atol

# no impact at maturity
Dν = Dual(Dual(ν₀,0,1), Dual(1,0,0))
zerogreeks = u(T, Ds, Dν, heston, EuropeanPut(T,K))

@test zerogreeks.value.partials[1] ≈ 0 atol=atol
@test zerogreeks.value.partials[2] ≈ 0 atol=atol
@test zerogreeks.partials[1].partials[1] ≈ 0 atol=atol
@test zerogreeks.partials[1].partials[2] ≈ 0 atol=atol
