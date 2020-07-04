include("../base.jl")
include("../params.jl")

using Plots
using LaTeXStrings

import PyPlot
pyplot()

t = 0.2
T = 0.5
K = 100
K₁ = 95
K₂ = 110
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

function diffop_plot(t, x, y, contract::AbstractContract)
    f(x,y) = diffop(t, x, y, heston, contract)
    surface(x, y, f)
    xlabel!(L"S_t")
    ylabel!(L"\nu_t")
    # savefig("DiffusionOperator"*string(typeof(contract))*".pdf")
end

nsteps = 51

x = range(80, stop=120, length=nsteps)
x_pow = range(8, stop=12, length=nsteps)
y = range(0.001, stop=0.025, length=nsteps)

theme(
    :default,
    colorbar=false,
    size=(800, 400),
    color=:ice,
)

d1 = diffop_plot(t, x, y, EuropeanPut(T,K))
d2 = diffop_plot(t, x, y, Strangle(T,K₁,K₂))
d3 = diffop_plot(t, x, y, PowerBinaryPut(T,K,n))
d4 = diffop_plot(t, x_pow, y, PowerCall(T,K,n))
