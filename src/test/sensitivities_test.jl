# test automatic differentiation against analytical solutions to undiscounted greeks / sensitivities
const BLACKSCHOLES_DELTA_PUT = -0.397681908481585
const BLACKSCHOLES_VEGA_PUT = 39.6952547477012
const BLACKSCHOLES_VANNA_PUT = -0.198476273738506
const BLACKSCHOLES_VOMMA_PUT = 5.95428821215518
const GENERALIZEDBLACKSCHOLES_DELTA_PUT = -0.40114740442307
const GENERALIZEDBLACKSCHOLES_VEGA_PUT = 41.463541599198
const GENERALIZEDBLACKSCHOLES_VANNA_PUT = 0.0799704029843022
const GENERALIZEDBLACKSCHOLES_VOMMA_PUT = -122.130268204048

# black-scholes
zs = Dual(Dual(s₀,1,0), Dual(0,0,0))
zσ = Dual(Dual(σ,0,1), Dual(1,0,0))
bssensitivities = u(0., zs, zσ, r, BlackScholes(), EuropeanPut(T,K))

@test bssensitivities.value.partials[1] ≈ BLACKSCHOLES_DELTA_PUT atol=atol
@test bssensitivities.value.partials[2] ≈ BLACKSCHOLES_VEGA_PUT atol=atol
@test bssensitivities.partials[1].partials[1] ≈ BLACKSCHOLES_VANNA_PUT atol=atol
@test bssensitivities.partials[1].partials[2] ≈ BLACKSCHOLES_VOMMA_PUT atol=atol

# generalized black-scholes
zν = Dual(Dual(ν₀,0,1), Dual(1,0,0))
gbssensitivities = u(0., zs, zν, heston, EuropeanPut(T,K))

@test gbssensitivities.value.partials[1] ≈ GENERALIZEDBLACKSCHOLES_DELTA_PUT atol=atol
@test gbssensitivities.value.partials[2] ≈ GENERALIZEDBLACKSCHOLES_VEGA_PUT atol=atol
@test gbssensitivities.partials[1].partials[1] ≈ GENERALIZEDBLACKSCHOLES_VANNA_PUT atol=atol
@test gbssensitivities.partials[1].partials[2] ≈ GENERALIZEDBLACKSCHOLES_VOMMA_PUT atol=atol

# no impact at maturity
zν = Dual(Dual(ν₀,0,1), Dual(1,0,0))
zerosensitivities = u(T, zs, zν, heston, EuropeanPut(T,K))

@test zerosensitivities.value.partials[1] ≈ 0 atol=atol
@test zerosensitivities.value.partials[2] ≈ 0 atol=atol
@test zerosensitivities.partials[1].partials[1] ≈ 0 atol=atol
@test zerosensitivities.partials[1].partials[2] ≈ 0 atol=atol
