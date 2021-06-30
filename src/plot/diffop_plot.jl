include("../base.jl")
include("../params.jl")

using Plots
using LaTeXStrings

import PyPlot
pyplot()

t = 0.7
T = 1
K = 100
H = 90
barrier_nsteps = 200

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

y = range(0.001, stop=0.025, length=nsteps)

theme(
    :default,
    colorbar=false,
    size=(800, 400),
    color=:ice,
)

put_x = range(80, stop=120, length=nsteps)
diffop_plot(t, x, y, EuropeanPut(T,K), DefaultState())

barrier_x = range(90, stop=120, length=nsteps)
diffop_plot(t, barrier_x, y, DownOutBarrierCall(T,K,H,barrier_nsteps), BarrierBreach(false))

lookback_x = range(100, stop=120, length=nsteps)
diffop_plot(t, lookback_x, y, FloatingLookbackPut(T), RunningMax(120))
