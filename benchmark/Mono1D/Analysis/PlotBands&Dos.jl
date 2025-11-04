include("../Mono1DParameters.jl")
using PyPlot
using Colors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))

Energies = range(-2.1, 2.1, 201)
RefDos = Ref.(Energies, 1)
thepath = range(-0.5, 0.5, 1501)

gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
soft_olive = "#a6a631"
clear_gray = "#c0c0c0"

bands = reduce(hcat, map(k -> H(k), thepath))
fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], color = gold, linewidth = 3)
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 751, 1501])
ax[1].set_xticklabels([-0.5, 0.0, 0.5])
ax[1].hlines(0.0, 1, 1501, color = :g, linewidth = 3, linestyle = :dashdot)
ax[1].hlines(1.5, 1, 1501, color = :b, linewidth = 3, linestyle = :dotted)
ax[1].hlines(1.99, 1, 1501, color = :r, linewidth = 3, linestyle = :dashed)
#ax[1].hlines([-7.72, 7.83, 10.05], 1, 1501, colors = [:k, :k, :k], linestyles = :dashed, linewidths = 2, label = "Crossing")
ax[1].set_ylim(-2.1, 2.1)
ax[1].set_xlim(0, 1501)
ax[1].set_ylabel("Energies (eV)")

ax[2].plot(RefDos, Energies, color = :black, linewidth = 3)
ax[2].fill_betweenx(Energies, RefDos, color = :black, alpha = 0.1)
ax[2].hlines(0.0, -1e-1, 2.5, color = :g, linewidth = 3, label = "Easy", linestyle = :dashdot)
ax[2].hlines(1.5, -1e-1, 2.5, color = :b, linewidth = 3, label = "Medium", linestyle = :dotted)
ax[2].hlines(1.99, -1e-1, 2.5, color = :r, linewidth = 3, label = "Hard", linestyle = :dashed)
#ax[2].hlines([-7.72, 7.83, 10.05], 0, 0.7, colors = [:k, :k, :k], linestyles = :dashed, linewidths = 2)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(-1e-2, 2.3)
ax[2].set_ylim(-2.1, 2.1)
#ax[2].set_ylim(-22.4, 23)
#ax[2].legend()

fig.legend(loc = "center right", bbox_to_anchor = (0.9, 0.3))
