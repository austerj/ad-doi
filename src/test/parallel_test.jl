# test that single/multi-threaded DOI estimates are identically distributed
using HypothesisTests: ApproximateTwoSampleKSTest, pvalue

# run estimators
single = [estimator(nsteps, npaths_run, heston, EuropeanPut(T,K), rng)[2] for i = 1:nruns]
parallel = [parallel_estimator(nsteps, npaths_run, heston, EuropeanPut(T,K))[2] for i = 1:nruns]

# use two-sample Kolmogorov-Smirnov test
@test pvalue(ApproximateTwoSampleKSTest(single, parallel)) > significance
