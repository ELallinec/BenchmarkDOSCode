include("../GrapheneParameters.jl")
using PyPlot
using Colors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
G = zeros(2)
M = 0.5 * [1, 1]
K = (1 / 3) * [2, 1]
tabpath = [G, M, K, G]
function segment(p1, p2, N)
	return [(1 - n / N) * p1 + n / N * p2 for n in 0:N]
end

function path(tabpath, N)
	return reduce(vcat, [(i == length(tabpath) - 1 ? segment(tabpath[i], tabpath[i+1], N) : segment(tabpath[i], tabpath[i+1], N)[1:end-1]) for i in eachindex(tabpath[1:end-1])])
end
Energies = range(-3.1, 3.1, 1001)
RefDos = dos_graphene_exact.(Energies)
thepath = path(tabpath, 500)

dark_navy = "#003f5c"
gold = "#ffa600"

bands = reduce(hcat, map(k -> real(eigen(H(k)).values), thepath))
fig, ax = subplots(1, 2, width_ratios = [2, 2])
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], color = gold, linewidth = 3)
ax[1].plot(bands[2, :], color = dark_navy, linewidth = 3)
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 501, 1001, 1501])
ax[1].set_xticklabels(["G", "M", "K", "G"])
ax[1].vlines([1, 501, 1001, 1501], -24, 24, colors = :black)
ax[1].hlines(2.0, 1, 1501, color = :g, linewidth = 3, linestyle = :dashdot)
ax[1].hlines(0.1, 1, 1501, color = :b, linewidth = 3, linestyle = :dotted)
ax[1].hlines(0.99, 1, 1501, color = :r, linewidth = 3, linestyle = :dashed)
#ax[1].hlines([-7.72, 7.83, 10.05], 1, 1501, colors = [:k, :k, :k], linestyles = :dashed, linewidths = 2, label = "Crossing")
ax[1].set_ylim(-3.1, 3.1)
ax[1].set_xlim(0, 1501)
ax[1].set_ylabel("Energies (eV)")

ax[2].plot(RefDos, Energies, color = :black, linewidth = 3)
ax[2].fill_betweenx(Energies, RefDos, color = :black, alpha = 0.1)
ax[2].hlines(2.0, -1e-1, 2, color = :g, linewidth = 3, label = "Easy", linestyle = :dashdot)
ax[2].hlines(0.1, -1e-1, 2, color = :b, linewidth = 3, label = "Medium", linestyle = :dotted)
ax[2].hlines(0.99, -1e-1, 2, color = :r, linewidth = 3, label = "Hard", linestyle = :dashed)
#ax[2].hlines([-7.72, 7.83, 10.05], 0, 0.7, colors = [:k, :k, :k], linestyles = :dashed, linewidths = 2)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(-1e-2, 1.3)
ax[2].set_ylim(-3.1, 3.1)
#ax[2].set_ylim(-22.4, 23)
#ax[2].legend()

fig.legend(loc = "center right", bbox_to_anchor = (0.9, 0.22))