using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
using JLD2

@load "benchmark/Mono1D/Results/ValuesBCDNoSmearing_N3:1500.jld2"
@load "benchmark/Mono1D/Results/ValuesPTR_N3:1500.jld2"
@load "benchmark/Mono1D/Results/ValuesLT_N3:1500.jld2"
@load "benchmark/Mono1D/Results/AllTimes.jld2"

for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(BCDtimes[i], abs.(BCDvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.scatter(LTtimes[i], abs.(LTvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT")
	ax.scatter(PTRtimes[i], abs.(PTRvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.scatter(IAItimes[i], abs.(allerrorsIAI[i]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI")
	yscale(:log)
	xlim(tabN[1], tabN[end])
	xlabel("Time (s)")
	ylabel("Relative error")
	ax.legend()
	savefig("benchmark/Mono1D/fig/Mono1D_RelErrorsVsTime_E=$(Energies[i]).png")
end