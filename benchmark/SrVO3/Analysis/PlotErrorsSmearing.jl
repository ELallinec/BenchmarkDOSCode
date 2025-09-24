using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 33)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2
path = "benchmark/SrVO3/ResultsSmearing/"


exDOSη, Energies, ηlist = load(path * "exDOS_Smearing_Multiple_final.jld2", "exDOSη", "Energies", "ηlist")
PTRvalues, PTRtimes = load("benchmark/SrVO3/Results/ValuesPTR_Multipleη_1e-5.jld2", "PTRval", "PTRtim")
BCDvalues, BCDtimes = load(path * "ValuesBCD_Multipleη_N3:200_final.jld2", "BCDvalues", "BCDtimes")
IAIvalues, IAItimes, IAINs = load("benchmark/SrVO3/Results/ValuesIAI_Multipleη_1e-5.jld2", "IAIval", "IAItim", "IAIN")

function premier_indice_complet(indices::Vector{Int}, indice_max::Int)
	sorted_indices = sort(unique(indices))  # On supprime les doublons au cas où
	n = length(sorted_indices)

	for i in 1:n
		start = sorted_indices[i]
		if start > indice_max
			continue  # inutile de vérifier des indices au-delà de la borne max
		end
		attendu = start:indice_max
		if all(x -> x in indices, attendu)
			return start
		end
	end

	return Inf# Aucun indice ne remplit la condition
end


tabindBCD = zeros(Int64, length(Energies), length(ηlist))
for ie in eachindex(Energies)
	for iη in eachindex(ηlist)
		for iN in 1:(size(BCDvalues, 3)-1)
			if abs(BCDvalues[ie, iη, iN+1] - BCDvalues[ie, iη, iN]) <= 1e-5
				println(BCDvalues[ie, iη, iN+1])
				println("E=$(Energies[ie]) ", "η=$(ηlist[iη]) ", "N= $(iN+2)")
				tabindBCD[ie, iη] = Int64(iN + 1)
				break
			end
		end
	end
end
#=
threshold = 5e-6
tempPTR = [[findall(k -> abs.(PTRvalues[i, j, k] - exDOSη[i, j]) <= threshold, 1:398) for j in eachindex(ηlist)] for i in eachindex(Energies)]
temp2PTR = [premier_indice_complet.(tempPTR[i], 398) for i in eachindex(Energies)]
PTRtime = zeros(length(Energies), length(ηlist))
for i in eachindex(Energies)
	for j in eachindex(ηlist)
		if isinf(temp2PTR[i][j])
			PTRtime[i, j] = Inf
		else
			PTRtime[i, j] = minimum(PTRtimes[i, j, temp2PTR[i][j]:end])
		end
	end
end


tempIAI = [[findall(k -> abs.(IAIvalues[i, j, k] - exDOSη[i, j]) <= threshold, eachindex(tolerances)) for j in eachindex(ηlist)] for i in eachindex(Energies)]
temp2IAI = [premier_indice_complet.(tempIAI[i], length(tolerances)) for i in eachindex(Energies)]
IAItime = zeros(length(Energies), length(ηlist))
for i in eachindex(Energies)
	for j in eachindex(ηlist)
		if isinf(temp2IAI[i][j])
			IAItime[i, j] = Inf
		else
			IAItime[i, j] = minimum(IAItimes[i, j, temp2IAI[i][j]:end])
		end
	end
end


tempBCD = [[findall(k -> abs.(BCDvalues[i, j, k] - exDOSη[i, j]) <= threshold, 1:198) for j in eachindex(ηlist)] for i in eachindex(Energies)]
temp2BCD = [premier_indice_complet.(tempBCD[i], 198) for i in eachindex(Energies)]
BCDtime = zeros(length(Energies), length(ηlist))
for i in eachindex(Energies)
	for j in eachindex(ηlist)
		if isinf(temp2BCD[i][j])
			BCDtime[i, j] = Inf
		else
			BCDtime[i, j] = minimum(BCDtimes[i, j, temp2BCD[i][j]:end])
		end
	end
end
=#
BCDtime = reduce(hcat, [[BCDtimes[ie, iη, tabindBCD[ie, iη]] for iη in eachindex(ηlist)] for ie in eachindex(Energies)])'
for i in eachindex(Energies)
	fig, ax = subplots()
	ax.plot(ηlist, PTRtimes[i, :], linestyle = :dashdot, color = (0, 158 / 255, 115 / 255), label = "PTR", linewidth = 3)
	ax.scatter(ηlist, PTRtimes[i, :], marker = "2", color = (0, 158 / 255, 115 / 255), s = 250)

	ax.plot(ηlist, IAItimes[i, :], "--k", label = "IAI", linewidth = 3)
	ax.scatter(ηlist, IAItimes[i, :], marker = "^", color = :black, s = 100)

	ax.plot(ηlist, BCDtime[i, :], linestyle = "-", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", linewidth = 3)
	ax.scatter(ηlist, BCDtime[i, :], marker = "v", color = (204 / 255, 121 / 255, 167 / 255), s = 100)

	#ax.plot(ηlist[9:14], log.(1 ./ ηlist[9:14]) .^ 3 / 10, linestyle = :dashed, color = :red, label = "Log3", linewidth = 3)
	#ax.plot(ηlist[3:7], ηlist[3:7] .^ (-3) / 1e4, linestyle = :dotted, color = :red, label = "1/η", linewidth = 3)
	ax.set_xscale(:log)
	ax.set_yscale(:log)
	ax.set_xlabel("η")
	ax.set_ylabel("Time (s)")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	#savefig("SrVO3_Time_vs_Eta_E=$(Energies[i]).pdf")
end

close("all")

plot(vcat(range(10.6, 11, 10), Energies), vcat(zeros(10), RefDos), "-k", linewidth = 5)
xlabel("Energies (eV)")
ylabel("D(E)")
vlines(12.44, -0.2, 6.2, linestyle = :dashdot, color = :green, lw = 4, label = "Easy")
vlines(13.29, -0.2, 6.2, linestyle = :dotted, color = :blue, lw = 4, label = "Medium")
vlines(11.55, -0.2, 6.2, linestyle = :dashed, color = :red, lw = 4, label = "Hard")

#text(11.9, 5.5, "E=12.44")
#text(11.01, 5.5, "E=11.55")
#text(12.75, 5.5, "E=13.29")
ylim(-0.1, 6.1)
xlim(10.6, 14.0)
legend(fontsize = 22, loc = "upper left")
