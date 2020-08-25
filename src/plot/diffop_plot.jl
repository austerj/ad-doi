include("../base.jl")
include("../params.jl")

using Plots
using LaTeXStrings

import PyPlot
pyplot()

t = 0.7
T = 1
K = 100
K1 = 95
K2 = 110
n = 2
    
function remove_grid!()
    ax = Plots.PyPlot.gca()
    ax.xaxis.pane.fill = false
    ax.yaxis.pane.fill = false
    ax.zaxis.pane.fill = false
    ax.xaxis.pane.set_edgecolor("w")
    ax.yaxis.pane.set_edgecolor("w")
    ax.zaxis.pane.set_edgecolor("w")
    ax.grid(false)
end

function diffop_plot(t, x, y, contract::AbstractContract, path_state::AbstractState)
    f(x,y) = diffop(t, x, y, heston, contract, path_state)
    surface(x, y, f)
    xlabel!(L"S_t")
    ylabel!(L"\nu_t")
    # savefig("DiffusionOperator"*string(typeof(contract))*".pdf")
end
nsteps = 51

x = range(80, stop=120, length=nsteps)
y = range(0.001, stop=0.025, length=nsteps)

theme(
    :default,
    colorbar=false,
    size=(800, 400),
    color=:ice,
)

d1 = diffop_plot(t, x, y, EuropeanPut(T,K), DefaultState())
d2 = diffop_plot(t, x, y, Strangle(T,K1,K2), DefaultState())
d3 = diffop_plot(t, x, y, PowerBinaryPut(T,K,n), DefaultState())
d4 = diffop_plot(t, x, y, PowerCall(T,K^n,n), DefaultState())
d5 = diffop_plot(t, x, y, FloatingLookbackPut(T), RunningMax(120))
