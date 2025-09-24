include("../SrVO3Parameters.jl")

using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (16, 9))
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

thepath = path(tabpath, 100)

bands = reduce(hcat, map(k -> real(eigen(H(k)).values), thepath))

fig, ax = subplots(1, 2, width_ratios = [3, 2])
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands', linewidth = 2)
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 101, 201, 301, 401])
ax[1].set_xticklabels(["Γ", "M", "X", "G", "R"])
ax[1].set_yticks([11, 12, 13, 14])
ax[1].vlines([0, 101, 201, 301, 401], 11, 14, colors = :black)
ax[1].set_ylim(11, 14)
ax[1].set_xlim(0, 401)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(RefDos, Energies, color = :black, linewidth = 1)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(0, 6)
ax[2].set_ylim(11, 14)
ax[2].fill_betweenx(Energies, RefDos, color = :black, alpha = 0.1)

# Draw hlines (now only one call is needed, but still u	se xmin/xmax for each subplot)
ax[1].hlines([12.40, 12.46], xmin = 0, xmax = 400, color = "red", linestyle = "--")
ax[2].hlines([12.40, 12.46], xmin = ax[2].get_xlim()[1], xmax = ax[2].get_xlim()[2], color = "red", linestyle = "--")



fig, ax = subplots()
ax.plot(bands', linewidth = 3)
ax.set_xlabel("k")
ax.set_xticks([1, 101, 201, 301, 401])
ax.set_xticklabels(["Γ", "M", "X", "G", "R"])
ax.set_yticks([11, 12, 13, 14])
ax.vlines([0, 101, 201, 301, 401], 11, 14, colors = :black)
ax.set_ylim(11, 14)
ax.set_xlim(0, 401)
ax.set_ylabel("Energies (eV)")