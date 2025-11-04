include("../GrapheneParameters.jl")
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))


tabk = range(-0.5, 0.5, 100)
tabk = Iterators.product(range(-0.5, 0.5, 100), range(-0.5, 0.5, 100))
tabk = [[k[1], k[2]] for k in tabk]
tabe1 = getindex.(eigvals.(H.(tabk)), 1)
tabe2 = getindex.(eigvals.(H.(tabk)), 2)

plot3D(getindex.(tabk, 1), getindex.(tabk, 2), tabe1)

PTRerror = [minimum(abs.(PTRvalues[i, :, iN] .- exDOS[i])) for i in eachindex(Energies), iN in 1:398]
PTRtime = [minimum(PTRtimes[i, :, iN]) for i in eachindex(Energies), iN in 1:398]
for i in eachindex(Energies)[2:end]
	fig, ax = subplots()
	ax.scatter(PTRtime[i, :], abs.(PTRerror[i, :]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR", s = 150)
	ax.scatter(IAItimes[i, :, :], abs.(IAIvalues[i, :, :] .- exDOS[i]) ./ (abs(exDOS[i]) == 0 ? 1 : abs(exDOS[i])), marker = "^", color = :black, label = "IAI", s = 100)
	ax.scatter(LTtimes[i, :], abs.(LTvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT", s = 150)
	ax.scatter(BCDtimes[i, :], abs.(BCDvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", s = 100)
	yscale(:log)
	xscale(:log)
	ylabel("Relative error")
	xlabel("Time (s)")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	savefig("benchmark/Graphene/fig/Graphene_ErrorsVsTime_E=$(Energies[i]).pdf")
end

close("all")
