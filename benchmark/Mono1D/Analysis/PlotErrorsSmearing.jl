using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2

@load "benchmark/Mono1D/Results/ValuesBCD_Smearing_N3:3165.jld2"
@load "benchmark/Mono1D/Results/ValuesPTR_Smearing_N3:3165.jld2"
@load "benchmark/Mono1D/Results/ErrorsIAI_Smearing.jld2"

for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(tabN, abs.(BCDvalues[:, i] .- exDOSη[i]) ./ (abs(exDOSη[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.scatter(tabN, abs.(PTRvalues[:, i] .- exDOSη[i]) ./ (abs(exDOSη[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.scatter(tabNIAI[i], abs.(allerrorsIAI[i]) ./ (abs(exDOSη[i])), marker = "^", color = :black, label = "IAI")
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("N")
	ylabel("Relative error")
	ax.legend()
	ax.set_title("$(Energies[i])")
	#savefig("benchmark/Mono1D/fig/Mono1D_RelErrorsVsN_SmearedDOS_E=$(Energies[i]).png")
end
