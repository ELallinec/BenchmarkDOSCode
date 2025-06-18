include("../GrapheneParameters.jl")
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 16)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))

x = y = range(-0.5, 0.49, 30)
BZ = DensityOfStates.BZ(d, 30)

tabh = (H.(BZ) + adjoint.(H.(BZ))) / 2
tabε = map(h -> eigvals(h)[2], tabh)
tabd1h = map(k -> DH(k)[2], BZ)

tempE = 2

temp = findall(i -> abs(tabε[i] .- tempE) <= 5e-4, CartesianIndices(tabε))
xs = x[[t[1] for t in temp]]
ys = y[[t[2] for t in temp]]
cs = contour(x, y, tabε, levels = [tempE], color = :b)
p = cs.collections[1].get_paths()[1]
v = p.vertices
xs = v[:, 1]
ys = v[:, 2]
temp1 = 10
for i in eachindex(x)
	for j in eachindex(y)
		h, d1h, d2h = DH([x[i], y[j]])
		temp = norm(-7e-3 * real.([tr(d1h[m] * exp(-((tabh[i, j] - tempE * I) / 0.5)^2)) for m in eachindex(d1h)]))
		if temp1 > temp && temp > 0
			temp1 = temp
		end
	end
end

using PyCall
for i in eachindex(x)
	for j in eachindex(y)
		h, d1h, d2h = DH([x[i], y[j]])
		temp = -7e-3 * real.([tr(d1h[m] * exp(-((tabh[i, j] - tempE * I) / 0.5)^2)) for m in eachindex(d1h)])

		arrow(x[i], y[j], temp[1], temp[2])
	end
end
xlabel(L"k_x")
ylabel(L"k_y")

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
fig,ax = subplots()
fig = figure()
ax = fig.add_subplot(111, projection="3d") # Correct way in Julia

plot_surface(tabx, taby, εm, linewidth = 0, cmap = :viridis)
plot_surface(tabx, taby, εp, linewidth = 0, cmap = :viridis)

