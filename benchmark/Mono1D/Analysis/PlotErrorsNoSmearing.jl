using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2

@load "benchmark/Mono1D/Results/ValuesBCD_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesPTR_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesLT_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesIAI_final.jld2"
PTRtime = [minimum(PTRtimes[i, :, iN]) for i in 1:5, iN in 1:1498]
PTRerror = [minimum(abs.(PTRvalues[i, :, iN] .- exDOS[i])) for i in 1:5, iN in 1:1498]
for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(3:1500, PTRerror[i, :] ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR", s = 150)
	ax.scatter(tabNIAI[i], abs.(allerrorsIAI[i]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI", s = 150)
	ax.scatter(3:1500, abs.(LTvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT", s = 150)
	ax.scatter(3:1500, abs.(BCDvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", s = 150)
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("N")
	ylabel("Relative error")
	ax.set_title("E=$(Energies[i])")
	if i != 5
		ax.legend(bbox_to_anchor = (0.2, 0, 0.4, 0.5))
	else
		ax.legend()
	end
	#savefig("benchmark/Mono1D/fig/Mono1D_RelErrorsVsN_E=$(Energies[i]).pdf")
end

tabind=[[argmin(abs.(BCDvalues[ie,:,iN].-exDOS[ie])) for iN in 1:1498] for ie in eachindex(Energies)]
temp=[[abs.(BCDvalues[ie,tabind[ie][iN],iN].-exDOS[ie]) for iN in 1:1498] for ie in eachindex(Energies)]
BCDtime = [[BCDtimes[ie,tabind[ie][iN],iN] for iN in 1:1498] for ie in eachindex(Energies)]
close("all")	
Energies = -2.5:0.008:2.5
plot(Energies, Ref.(Energies, 1), "-k", linewidth = 5)
xlabel("Energies")
ylabel("D(E)")
vlines(0.0, -0.2, 3.0, linestyle = :dashdot, color = :green, lw = 4, label = "Easy")
vlines(1.5, -0.2, 3, linestyle = :dotted, color = :blue, lw = 4, label = "Medium")
vlines(1.99, -0.2, 3, linestyle = :dashed, color = :red, lw = 4, label = "Hard")
xlim(-2.2, 2.2)
ylim(-0.1, 2.6)
legend(fontsize = 28, loc = "best", bbox_to_anchor = (-0.07, 0.5, 0.5, 0.5))
