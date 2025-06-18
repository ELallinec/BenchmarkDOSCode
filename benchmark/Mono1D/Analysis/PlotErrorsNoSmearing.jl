using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
using JLD2

@load "benchmark/Mono1D/Results/AllErrors.jld2"

for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(3:1500, abs.(BCDvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.scatter(3:1500, abs.(LTvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT")
	ax.scatter(3:1500, abs.(PTRvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.scatter(tabNIAI[i], abs.(allerrorsIAI[i]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI")
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("N")
	ylabel("Relative error")
	ax.legend()
	savefig("benchmark/Mono1D/fig/Mono1D_RelErrorsVsN_E=$(Energies[i]).png")
end