include("../GrapheneParameters.jl")
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))

x = y = range(-0.5, 0.49, 51)
BZ = [[i[1], i[2]] for i in Iterators.product(x, y)]

tabh = (tb_graphene.(BZ) + adjoint.(tb_graphene.(BZ))) / 2
tabε = map(h -> eigvals(h)[2], tabh)
tempE = 2

levels = [tempE]
cs = contour(x, y, tabε, levels = levels, colors = ["red"], linewidths = [5])
xs = ys = range(-0.5, 0.49, 17)
smallBZ = [[i[1], i[2]] for i in Iterators.product(xs, ys)]
using PyCall
for ind in smallBZ[1:end-1]
	h, d1h, d2h = DH([ind[1], ind[2]])
	temp = -7e-3 * real.([tr(d1h[m] * exp(-((h - tempE * I) / 0.5)^2)) for m in eachindex(d1h)])

	arrow(ind[1], ind[2], temp[1], temp[2], head_width = 0.35 * norm(temp), width = 0.003, color = :k)
end
plot([], [], color = "red", lw = 5, label = "Fermi Level")  # empty plot with color + label
plot([], [], c = "black", marker = "\$\\rightarrow\$", linestyle = "none", markersize = 40, label = "Deformation")
xlabel(L"k_x")
ylabel(L"k_y")
legend(loc = "upper left", fontsize = 22)

include("../GrapheneParameters.jl")

tabx = range(-0.5, 0.5, 100)
taby = range(-0.5, 0.5, 100)
meshing = [[k[1], k[2]] for k in Iterators.product(tabx, taby)]
εm = getindex.(eigvals.(H.(meshing)), 1)
εp = getindex.(eigvals.(H.(meshing)), 2)




plot_surface(tabx, taby, εp, linewidth = 0, cmap = :viridis)
locator_params(nbins = 3)
xlabel(L"k_x")
ylabel(L"k_y")
zlabel(latexstring("\\varepsilon_{\\pm}"))
fig, ax = subplots()
fig = figure()
ax = fig.add_subplot(111, projection = "3d") # Correct way in Julia

plot_surface(tabx, taby, εm, linewidth = 0, cmap = :viridis)
plot_surface(tabx, taby, εp, linewidth = 0, cmap = :viridis)

