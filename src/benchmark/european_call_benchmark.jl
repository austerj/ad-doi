include("../base.jl")
include("../params.jl")
include("benchmarks.jl")

using Printf

# variations based on parameter set A
ν₀_A = [0.01, 0.04, 0.16]
K_A = [90, 95, 100, 105, 110]
T = 1
nsteps = 200

# number of paths per estimates and number of estimates
npaths = 10000
nestimates = 50

function print_table()
    println("\\begin{tabular}{ccrccrcc}")
    println("\\toprule")
    print("&&\\multicolumn{2}{c}{MC} & \\,&\\multicolumn{2}{c}{AD DOI}")
    println(" \\\\")
    println("\\cmidrule{3-4} \\cmidrule{6-7}")
    println("K & \$\\nu_0\$ & Mean & Conf.Interval && Mean & Conf.Interval \\\\")
    for j in 1:size(K_A)[1]
        println("\\midrule")
        K = K_A[j]
        @printf("\\multirow{3}{*}{%d}\n", K)
        for i in 1:size(ν₀_A)[1]
            heston_A = Heston(s₀, ν₀_A[i], r, κ, θ, ξ, ρ)
            contract = EuropeanCall(T,K)
            mc_mean, mc_conf, doi_mean, doi_conf = confidence_benchmark(nsteps, npaths, heston_A, contract, nestimates)

            @printf(" & %.2f", ν₀_A[i])
            @printf(" & %.4f & (%.4f,%.4f)", mc_mean[1], mc_conf[1][1], mc_conf[2][1])
            @printf(" && %.4f & (%.4f,%.4f)", doi_mean[1], doi_conf[1][1], doi_conf[2][1])
            println(" \\\\")
        end
    end
    println("\\bottomrule")
    println("\\end{tabular}")
end

print_table()

