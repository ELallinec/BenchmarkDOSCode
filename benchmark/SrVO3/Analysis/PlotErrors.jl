using Pkg: Pkg
push!(LOAD_PATH, "src")


#using DensityOfStates
using JLD2
using PyPlot
using PlotErrors
include("../SrVO3Parameters.jl")
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
#PyPlot.tight_layout()

@load "benchmark/SrVO3/Results/Alltime+errors.jld2"
ind = findall(i -> i <= tabN[end], tabNIAI)[end]
inter = tabNIAI[1:ind]
#savepath="benchmark/Graphene/fig/nvserrorsgraphene.pdf"
savepath = false
fig, ax =
	PlotErrors.scatterresults42(
		tabN,
		tabNIAI[1:ind],
		errorsIAIL1[1:ind] / norm(RefDos, 1),
		errorsLTL1 / norm(RefDos, 1),
		errorsBCDL1 / norm(RefDos, 1),
		errorsPTRL1 / norm(RefDos, 1),
		savepath,
		1:1:length(tabN),
		"N",
		latexstring("\$ \\ell^{1}\$ error"),
		"best",
		"linear",
		"log";
	)



fig, ax =
	PlotErrors.scatterresults42(
		tabN,
		tabNIAI[1:ind],
		errorsIAILinf[1:ind] / norm(RefDos, Inf),
		errorsLTLinf / norm(RefDos, Inf),
		errorsBCDLinf / norm(RefDos, Inf),
		errorsPTRLinf / norm(RefDos, Inf),
		savepath,
		1:1:length(tabN),
		"N",
		latexstring("\$ \\ell^{\\infty}\$ error"),
		"best",
		"linear",
		"log";
	)


fig, ax = subplots()
ax.scatter(timeLinfIAI / length(Energies), errorsIAILinf / maximum(RefDos), label = "IAI", color = :black, s = 35, marker = :^)
ax.scatter(timeLT[1:1:length(timeLT)] / length(Energies), errorsLTLinf[1:1:100] / maximum(RefDos), label = "LT", color = (230 / 255, 159 / 255, 0), s = 100, marker = :+)
ax.scatter(timeBCD[1:1:length(timeBCD)] / length(Energies), (errorsBCDLinf[eachindex(timeBCD)]) / maximum(RefDos), label = "BCD", color = (204 / 255, 121 / 255, 167 / 255), s = 100, marker = "2")
ax.scatter(timePTR[1:1:length(timePTR)] / length(Energies), errorsPTRLinf[1:1:100] / maximum(RefDos), label = "PTR", color = (0, 158 / 255, 115 / 255), s = 35, marker = :v)
yscale(:log)
xscale(:log)
ylim(5e-4, 5)

xlabel("Time (s)")
ylabel(latexstring("\$ \\ell^{\\infty}\$ error"))
legend()

fig, ax = subplots()

ax.scatter(timeL1IAI / length(Energies), errorsIAIL1 / norm(RefDos, 1), label = "IAI", color = :black, s = 35, marker = :^)
ax.scatter(timeLT[1:1:length(timeLT)] / length(Energies), errorsLTL1[1:1:100] / norm(RefDos, 1), label = "LT", color = (230 / 255, 159 / 255, 0), s = 100, marker = :+)
ax.scatter(timeBCD[1:1:length(timeBCD)] / length(Energies), (errorsBCDL1[eachindex(timeBCD)]) / norm(RefDos, 1), label = "BCD", color = (204 / 255, 121 / 255, 167 / 255), s = 100, marker = "2")
ax.scatter(timePTR[1:1:length(timePTR)] * 2 / length(Energies), errorsPTRL1[1:1:100] / norm(RefDos, 1), label = "PTR", color = (0, 158 / 255, 115 / 255), s = 35, marker = :v)

yscale(:log)
xscale(:log)
ylim(1e-4, 5)

xlabel("Time (s)")
ylabel(latexstring("\$ \\ell^{1}\$ error"))
legend()