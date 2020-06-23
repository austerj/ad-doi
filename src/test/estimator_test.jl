# test estimator against numerical integration (optByHestonNI) in MATLAB
const NUMERICALINTEGRATION_PUT = 12.0190225503127
const NUMERICALINTEGRATION_CALL = 15.9400786350804

# relative error tolerance of 0.1%
err_tol = 1e-3

# mean relative error ~1e-5 with nsteps=50, npaths=5000 based on 1000 samples for call
@test estimator(nsteps, npaths, heston, EuropeanPut(T,K), rng)[2] ≈ NUMERICALINTEGRATION_PUT atol=err_tol*NUMERICALINTEGRATION_PUT
@test estimator(nsteps, npaths, heston, EuropeanCall(T,K), rng)[2] ≈ NUMERICALINTEGRATION_CALL atol=err_tol*NUMERICALINTEGRATION_CALL
