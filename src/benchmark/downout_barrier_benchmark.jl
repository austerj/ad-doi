include("../base.jl")
include("benchmarks.jl")

using Printf

# parameter sets B1-B3
heston_B = Heston(100, .04, .04, .6, .04, .2, -.8)
T_B = 0.5
K_B = [90, 100, 110]
H_B = [92, 95, 98]
nsteps = Int(T_B/0.004)

# number of paths per estimates and number of estimates
npaths = 10000
nestimates = 50

# statistics from Coskun-Korn (2018), [mean, left-confidence, right-confidence]
coskun_korn = [
    # K = 90
    [
        [10.7489, 10.7466, 10.7512],  # H = 92
        [8.3899, 8.3847, 8.3951],  # H = 95
        [4.7526, 4.7454, 4.7598],  # H = 98
    ],
    # K = 100
    [
        [5.6910, 5.6849, 5.6971],  # H = 92
        [4.6761, 4.6724, 4.6799],  # H = 95
        [2.8158, 2.8147, 2.8168],  # H = 98
    ],
    # K = 110
    [
        [2.0608, 2.0538, 2.0679],  # H = 92
        [1.7750, 1.7678, 1.7823],  # H = 95
        [1.1523, 1.1459, 1.1587],  # H = 98
    ]
]

function print_table()
    println("\\begin{tabular}{ccrccrccrcc}")
    println("\\toprule")
    print("&&\\multicolumn{2}{c}{MC} & \\,&\\multicolumn{2}{c}{AD DOI} & \\,&\\multicolumn{2}{c}{Coskun-Korn}")
    println(" \\\\")
    println("\\cmidrule{3-4} \\cmidrule{6-7} \\cmidrule{9-10}")
    println("K & H & Mean & Conf.Interval && Mean & Conf.Interval && Mean & Conf.Interval \\\\")
    for j in 1:3
        println("\\midrule")
        K = K_B[j]
        @printf("\\multirow{3}{*}{%d}\n", K)
        for i in 1:3
            H = H_B[i]
            contract = DownOutBarrierCall(T_B,K,H,nsteps)
            mc_mean, mc_conf, doi_mean, doi_conf = confidence_benchmark(nsteps, npaths, heston_B, contract, nestimates)

            @printf(" & %d", H)
            @printf(" & %.4f & (%.4f,%.4f)", mc_mean[1], mc_conf[1][1], mc_conf[2][1])
            @printf(" && %.4f & (%.4f,%.4f)", doi_mean[1], doi_conf[1][1], doi_conf[2][1])
            @printf(" && %.4f & (%.4f,%.4f)", coskun_korn[j][i]...)
            println(" \\\\")
        end
    end
    println("\\bottomrule")
    println("\\end{tabular}")
end

print_table()

