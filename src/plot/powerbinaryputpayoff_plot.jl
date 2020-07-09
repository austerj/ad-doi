include("../base.jl")
include("../params.jl")

using Plots
using ColorSchemes
using LaTeXStrings

pyplot()

x = range(0, stop=120, length=10000)

theme(
    :default,
    size=(800, 250),
    xminorticks=4,
    yminorticks=3,
    xgrid=false,
    ygrid=false,
    legend=false
)
colors = [ColorSchemes.ice[i] for i in 64:128:192]
alpha = 0.8

p = plot()
contract = PowerBinaryPut(1, 100, 2);
payoff = [h(s, contract) for s in x]
plot!(x, payoff, color=colors[1], alpha=alpha)
xlabel!(L"S_t")

savefig("PowerBinaryPutPayoff.pdf")
