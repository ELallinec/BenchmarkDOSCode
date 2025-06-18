using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
using JLD2
path = "benchmark/Graphene/Results"
@load "benchmark/Graphene/Results/ValuesIAI.jld2"
@load "benchmark/Graphene/Results/ErrorsIAI.jld2"
@load "benchmark/Graphene/Results/ValuesLT_N3:150.jld2"
@load "benchmark/Graphene/Results/ValuesBCDNoSmearing_N3:150.jld2"
@load "benchmark/Graphene/Results/ValuesPTR_N3:150.jld2"

for i in eachindex(Energies)[2:end]
	fig, ax = subplots()
	ax.scatter(tabN, abs.(BCDvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.scatter(tabN, abs.(LTvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT")
	ax.scatter(tabN, abs.(PTRvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.scatter(tabNIAI[i], abs.(allerrorsIAI[i]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI")
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("N")
	ylabel("Relative error")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	savefig("benchmark/Graphene/fig/Graphene_RelErrorsVsN_E=$(Energies[i]).png")
end