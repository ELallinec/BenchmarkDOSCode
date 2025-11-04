include("../SrVO3Parameters.jl")

using PyPlot
using Colors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
G = zeros(3)
M = 0.5 * [1, 1, 0]
X = 0.5 * [0, 1, 0]
R = 0.5 * [1, 1, 1]
tabpath = [G, M, X, G, R]
function segment(p1, p2, N)
	return [(1 - n / N) * p1 + n / N * p2 for n in 0:N]
end

function path(tabpath, N)
	return reduce(vcat, [(i == length(tabpath) - 1 ? segment(tabpath[i], tabpath[i+1], N) : segment(tabpath[i], tabpath[i+1], N)[1:end-1]) for i in eachindex(tabpath[1:end-1])])
end

thepath = path(tabpath, 500)

bands = reduce(hcat, map(k -> real(eigen(H(k)).values), thepath))
@load "benchmark/SrVO3/Results/RefDos.jld2"
gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
soft_olive = "#a6a631"
clear_gray = "#c0c0c0"
fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], color = gold, linewidth = 3)
ax[1].plot(bands[2, :], color = dark_navy, linewidth = 3)
ax[1].plot(bands[3, :], color = rust_orange, linewidth = 3)
ax[1].hlines(12.44, 0, 2001, linewidth = 3, linestyle = :dashdot, color = :g, label = "Easy")
ax[1].hlines(13.29, 0, 2001, linewidth = 3, linestyle = :dotted, color = :b, label = "Medium")
ax[1].hlines(11.57, 0, 2001, linewidth = 3, linestyle = :dashed, color = :r, label = "Hard")
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 501, 1001, 1501, 2001])
ax[1].set_xticklabels(["Γ", "M", "X", "G", "R"])
ax[1].set_yticks([11, 12, 13, 14])
ax[1].vlines([1, 501, 1001, 1501, 2001], 11, 14, colors = :black)
ax[1].set_ylim(11.4, 13.9)
ax[1].set_xlim(1, 2001)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(RefDos, Energies, color = :black, linewidth = 3)
ax[2].hlines(12.44, 0, 2001, linewidth = 3, linestyle = :dashdot, color = :g)
ax[2].hlines(13.29, 0, 2001, linewidth = 3, linestyle = :dotted, color = :b)
ax[2].hlines(11.57, 0, 2001, linewidth = 3, linestyle = :dashed, color = :r)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(0, 6)
ax[2].set_ylim(11.4, 13.9)
ax[2].fill_betweenx(Energies, RefDos, color = :black, alpha = 0.1)
fig.legend(loc = "center right", bbox_to_anchor = (0.9, 0.3))

