module PlotErrors
using PyPlot


"""

	plotresults3(points, tabLT,tabBCD,tabSmearing,interval;kwargs...)

Plot 3 arrays on the same grid, often corresponding to the errors or times for BCD, Smearing and LT.
"""
function plotresults3(points, tabLT, tabBCD, tabSmearing, sav = false, interval = 1:length(points), xlabel = " ", ylabel = " ", loc = "best", xscale = "linear", yscale = "linear")

	fig = plot(points[interval], tabLT[interval], label = "LT", linestyle = :dashed, color = (230 / 255, 159 / 255, 0))
	scatter(points[interval], tabLT[interval], label = "", marker = :x, color = (230 / 255, 159 / 255, 0))

	plot(points[interval], tabBCD[interval], label = "BCD", linestyle = :solid, color = (204 / 255, 121 / 255, 167 / 255))
	scatter(points[interval], tabBCD[interval], label = "", marker = :x, color = (204 / 255, 121 / 255, 167 / 255))

	plot(points[interval], tabSmearing[interval], label = "PTR", linestyle = :dotted, color = (0, 158 / 255, 115 / 255))
	scatter(points[interval], tabSmearing[interval], label = "", marker = :x, color = (0, 158 / 255, 115 / 255))
	PyPlot.yscale(yscale)
	PyPlot.xscale(xscale)
	PyPlot.xlabel("$(xlabel)")
	PyPlot.ylabel("$(ylabel)")
	PyPlot.legend(loc = loc)
	if sav != false
		PyPlot.savefig(sav)
	end
	return fig
end

"""

	plotresults4(points,tabIAI,tabLT,tabBCD,tabSmearing,interval;kwargs...)

Plot 4 arrays on the same grid, often corresponding to the errors or times for each of our 4 methods.
"""
function plotresults4(points, tabIAI, tabLT, tabBCD, tabSmearing, sav = false, interval = 1:length(points), xlabel = " ", ylabel = " ", loc = "best", xscale = "linear", yscale = "linear")
	fig = plot(points[interval], tabIAI[interval], label = "IAI", linestyle = :dashdot, color = :black)
	scatter(points[interval], tabIAI[interval], label = "", marker = :x, color = :black)

	plot(points[interval], tabLT[interval], label = "LT", color = (230 / 255, 159 / 255, 0), linestyle = :dashed)
	scatter(points[interval], tabLT[interval], label = "", marker = :x, color = (230 / 255, 159 / 255, 0))

	plot(points[interval], tabBCD[interval], label = "BCD", linestyle = :solid, color = (204 / 255, 121 / 255, 167 / 255))
	scatter(points[interval], tabBCD[interval], label = "", marker = :x, color = (204 / 255, 121 / 255, 167 / 255))

	plot(points[interval], tabSmearing[interval], label = "PTR", linestyle = :dotted, color = (0, 158 / 255, 115 / 255))
	scatter(points[interval], tabSmearing[interval], label = "", marker = :x, color = (0, 158 / 255, 115 / 255))
	PyPlot.yscale(yscale)
	PyPlot.xscale(xscale)
	PyPlot.xlabel("$(xlabel)")
	PyPlot.ylabel("$(ylabel)")
	PyPlot.legend(loc = loc)
	if sav != false
		PyPlot.savefig(sav)
	end
	return fig
end


function plotresults42(points, pointsIAI, tabIAI, tabLT, tabBCD, tabSmearing, sav = false, interval = 1:length(points), xlabel = " ", ylabel = " ", loc = "best", xscale = "linear", yscale = "linear")
	fig, ax = subplots()
	ax.plot(pointsIAI, tabIAI, label = "IAI", linestyle = :dashdot, color = :black)
	ax.scatter(pointsIAI, tabIAI, label = "", marker = :x, color = :black)

	ax.plot(points[interval], tabLT[interval], label = "LT", color = (230 / 255, 159 / 255, 0), linestyle = :dashed)
	ax.scatter(points[interval], tabLT[interval], label = "", marker = :x, color = (230 / 255, 159 / 255, 0))

	ax.plot(points[interval], tabBCD[interval], label = "BCD", linestyle = :solid, color = (204 / 255, 121 / 255, 167 / 255))
	ax.scatter(points[interval], tabBCD[interval], label = "", marker = :x, color = (204 / 255, 121 / 255, 167 / 255))

	ax.plot(points[interval], tabSmearing[interval], label = "PTR", linestyle = :dotted, color = (0, 158 / 255, 115 / 255))
	ax.scatter(points[interval], tabSmearing[interval], label = "", marker = :x, color = (0, 158 / 255, 115 / 255))
	PyPlot.yscale(yscale)
	PyPlot.xscale(xscale)
	PyPlot.xlabel("$(xlabel)")
	PyPlot.ylabel("$(ylabel)")
	PyPlot.legend(loc = loc)
	if sav != false
		PyPlot.savefig(sav)
	end
	return fig, ax
end
function scatterresults42(points, pointsIAI, tabIAI, tabLT, tabBCD, tabSmearing, sav = false, interval = 1:length(points), xlabel = " ", ylabel = " ", loc = "best", xscale = "linear", yscale = "linear")
	fig, ax = subplots()
	#ax.plot(pointsIAI, tabIAI, label="IAI",linestyle=:dashdot,color=:black)
	ax.scatter(pointsIAI, tabIAI, label = "IAI", color = :black, s = 35, marker = :^)

	#ax.plot(points[interval],tabLT[interval], label="LT",color=(230/255,159/255,0),linestyle=:dashed)
	ax.scatter(points[interval], tabLT[interval], label = "LT", color = (230 / 255, 159 / 255, 0), s = 100, marker = :+)

	#ax.plot(points[interval],tabBCD[interval],label="BCD",linestyle=:solid,color=(204/255,121/255,167/255))
	ax.scatter(points[interval], tabBCD[interval], label = "BCD", color = (204 / 255, 121 / 255, 167 / 255), s = 100, marker = "2")

	#ax.plot(points[interval],tabSmearing[interval],label="PTR",linestyle=:dotted,color=(0,158/255,115/255))
	ax.scatter(points[interval], tabSmearing[interval], label = "PTR", color = (0, 158 / 255, 115 / 255), s = 35, marker = :v)
	PyPlot.yscale(yscale)
	PyPlot.xscale(xscale)
	PyPlot.xlabel("$(xlabel)")
	PyPlot.ylabel("$(ylabel)")
	PyPlot.legend(loc = loc)
	if sav != false
		PyPlot.savefig(sav)
	end
	return fig, ax
end



function plotthreshold4(points, tabIAI, tabLT, tabBCD, tabSmearing, sav = false, interval = 1:length(points), xlabel = " ", ylabel = " ", loc = "best", xscale = "linear", yscale = "linear")
	fig = plot(tabIAI[interval], points[interval], label = "IAI", linestyle = :dashdot, color = :black)
	scatter(tabIAI[interval], points[interval], label = "", marker = :x, color = :black)

	plot(tabLT[interval], points[interval], label = "LT", color = (230 / 255, 159 / 255, 0), linestyle = :dashed)
	scatter(tabLT[interval], points[interval], label = "", marker = :x, color = (230 / 255, 159 / 255, 0))

	plot(tabBCD[interval], points[interval], label = "BCD", linestyle = :solid, color = (204 / 255, 121 / 255, 167 / 255))
	scatter(tabBCD[interval], points[interval], label = "", marker = :x, color = (204 / 255, 121 / 255, 167 / 255))

	plot(tabSmearing[interval], points[interval], label = "PTR", linestyle = :dotted, color = (0, 158 / 255, 115 / 255))
	scatter(tabSmearing[interval], points[interval], label = "", marker = :x, color = (0, 158 / 255, 115 / 255))
	PyPlot.yscale(yscale)
	PyPlot.xscale(xscale)
	PyPlot.xlabel("$(xlabel)")
	PyPlot.ylabel("$(ylabel)")
	PyPlot.legend(loc = loc)
	if sav != false
		PyPlot.savefig(sav)
	end
	return fig
end

end
