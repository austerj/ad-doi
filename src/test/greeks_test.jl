# test automatic differentiation against analytical solutions to undiscounted greeks
const BLACKSCHOLES_DELTA_PUT = -0.212507002995933
const BLACKSCHOLES_VEGA_PUT = 29.5011840860567
const BLACKSCHOLES_VANNA_PUT = -0.924570912386131
const BLACKSCHOLES_VOMMA_PUT = 76.4437614171988
const GENERALIZEDBLACKSCHOLES_DELTA_PUT = -0.291290293982836
const GENERALIZEDBLACKSCHOLES_VEGA_PUT = 36.4943404584727
const GENERALIZEDBLACKSCHOLES_VANNA_PUT = -0.224847933750906
const GENERALIZEDBLACKSCHOLES_VOMMA_PUT = -91.6947601203105

# black-scholes greeks
Ds = Dual(Dual(s₀,1,0), Dual(0,0,0))
Dσ = Dual(Dual(σ,0,1), Dual(1,0,0))
state = BlackScholesState(0, Ds, Dσ, r)
bsgreeks = u(state, EuropeanPut(T,K))

@test bsgreeks.value.partials[1] ≈ BLACKSCHOLES_DELTA_PUT atol=atol
@test bsgreeks.value.partials[2] ≈ BLACKSCHOLES_VEGA_PUT atol=atol
@test bsgreeks.partials[1].partials[1] ≈ BLACKSCHOLES_VANNA_PUT atol=atol
@test bsgreeks.partials[1].partials[2] ≈ BLACKSCHOLES_VOMMA_PUT atol=atol

# generalized black-scholes
Dν = Dual(Dual(ν₀,0,1), Dual(1,0,0))
gbsgreeks = u(0, Ds, Dν, heston, EuropeanPut(T,K))

@test gbsgreeks.value.partials[1] ≈ GENERALIZEDBLACKSCHOLES_DELTA_PUT atol=atol
@test gbsgreeks.value.partials[2] ≈ GENERALIZEDBLACKSCHOLES_VEGA_PUT atol=atol
@test gbsgreeks.partials[1].partials[1] ≈ GENERALIZEDBLACKSCHOLES_VANNA_PUT atol=atol
@test gbsgreeks.partials[1].partials[2] ≈ GENERALIZEDBLACKSCHOLES_VOMMA_PUT atol=atol
