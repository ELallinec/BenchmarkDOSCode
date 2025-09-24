using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
using JLD2

@load "benchmark/Graphene/Results/ValuesIAI.jld2"

temp = [[abs.(IAIvalues[i, j, :] .- exDOSη[i, j]) for j in eachindex(ηlist)] for i in eachindex(Energies)]
for i in [2, 5, 6]
	fig, ax = subplots()
	for itol in eachindex(tolerances)
		ax.plot(Int64.(round.(IAINs[i, :, itol] .^ (1 / 2))), abs.(IAIvalues[i, :, itol] .- exDOS[i]) / exDOS[i], label = "tol=$(tolerances[itol])", "--", marker = :x)
	end
	ax.plot(tabNIAI[i], allerrorsIAI[i] / exDOS[i])
	ax.set_yscale(:log)
	ax.set_xlabel("N")
	ax.set_ylabel("Error")
	ax.set_title("E=$(Energies[i])")
	ax.legend(fontsize = 15, ncol = 3)
end


close("all")
