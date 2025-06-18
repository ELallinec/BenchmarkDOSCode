using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
using JLD2

@load "benchmark/Graphene/Results/ValuesIAI.jld2"

temp = [[abs.(IAIvalues[i, j, :] .- exDOSη[i, j]) for j in eachindex(ηlist)] for i in eachindex(Energies)]

temp2 = []
tempN2 = []
for i in eachindex(Energies)
	temp3 = []
	tempN3 = []
	for j in eachindex(ηlist)
		temp4 = findall(index -> temp[i][j][index] <= 1e-3, eachindex(tolerances))
		push!(temp3, temp4)
		if !isempty(temp4)
			push!(tempN3, minimum(IAINs[i, j, temp4]))
		else
			push!(tempN3, Inf)
		end

		#println(temp3)
	end
	push!(tempN2, tempN3)
	push!(temp2, temp3)
end


testtt = findall(index -> temp[1][1][index] <= 1e-3, eachindex(tolerances))
