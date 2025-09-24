using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2
path = "benchmark/Graphene/Results"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesIAI_final.jld2"
#@load "benchmark/Graphene/Results/ErrorsIAI.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesLT_N3:400_final.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesBCD_N3:400_final.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesPTR_N3:400_final.jld2"


PTRvalue = [[argmin(abs.(PTRvalues[ie, :, iN] .- exDOS[ie])) for iN in 1:398] for ie in eachindex(Energies)]
PTRerror = [[abs(PTRvalues[ie, PTRvalue[ie][iN], iN] - exDOS[ie]) for iN in 1:398] for ie in eachindex(Energies)]
PTRtime = [[PTRtimes[ie, PTRvalue[ie][iN], iN] for iN in 1:398] for ie in eachindex(Energies)]


LTvalue = [[argmin(abs.(LTvalues[ie, iN] .- exDOS[ie])) for iN in 1:398] for ie in eachindex(Energies)]
LTerror = [[abs(LTvalues[ie, LTvalue[ie][iN], iN] - exDOS[ie]) for iN in 1:398] for ie in eachindex(Energies)]
LTtime = [[LTtimes[ie, LTvalue[ie][iN], iN] for iN in 1:398] for ie in eachindex(Energies)]

for i in eachindex(Energies)[2:end]
	fig, ax = subplots()
	ax.scatter(tabN, abs.(PTRvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i]) == 0 ? 1 : abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR", s = 150)
	ax.scatter(tabNIAI[i], abs.(allerrorsIAI[i]) ./ (abs(exDOS[i]) == 0 ? 1 : abs(exDOS[i])), marker = "^", color = :black, label = "IAI", s = 100)
	ax.scatter(tabN, abs.(LTvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i]) == 0 ? 1 : abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT", s = 150)
	ax.scatter(tabN, abs.(BCDvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i]) == 0 ? 1 : abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", s = 100)
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("N")
	ylabel("Relative error")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	#savefig("benchmark/Graphene/fig/Graphene_RelErrorsVsN_E=$(Energies[i]).png")
end
close("all")

plot(-3.5:0.008:3.5, dos_graphene_exact.(-3.5:0.008:3.5), "-k", linewidth = 5)
xlabel("Energies")
ylabel("D(E)")
vlines(2.0, -0.2, 1.4, linestyle = :dashdot, color = :green, lw = 4, label = "Easy")
vlines(0.1, -0.2, 1.4, linestyle = :dotted, color = :blue, lw = 4, label = "Medium")
vlines(0.99, -0.2, 1.4, linestyle = :dashed, color = :red, lw = 4, label = "Hard")
#text(-0.65, 1.0, "E=0.1", fontsize = 25)
#text(0.67, 1.0, "E=0.99", fontsize = 25)
#text(2.05, 1.0, "E=2.0", fontsize = 25)
ylim(-0.05, 1.1)
xlim(-3.2, 3.2)
legend(fontsize = 28, loc = "upper left")
