using Pkg
Pkg.activate(".")
using JLD2
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (16, 9))

@load "benchmark/SrVO3/Results/ValuesIAI_Multipleη.jld2"

tempN = Int64.(round.(IAINs[2, :, :] .^ (1 / 3)))
utempN = unique(tempN)
indtempN = [findall(ind -> tempN[ind] == N, CartesianIndices(tempN)) for N in utempN]

temp1 = []
for i in eachindex(indtempN)
	indη = sort(unique(getindex.(indtempN[i], 1)))
	temp = []
	for j in eachindex(indη)
		push!(temp, findall(k -> indη[j] == k[1], indtempN[i]))
	end
	push!(temp1, temp)
end
tempvalues = IAIvalues[2, :, :]
tempN = Int64.(round.(IAINs[2, :, :] .^ (1 / 3)))
temp = []
for (icol, col) in enumerate(eachcol(tempN))
	utempN = unique(col)
	temp1 = []
	for N in utempN
		indicesN = findall(i -> col[i] == N, eachindex(col))
		push!(temp1, ηlist[argmin(tempvalues[indicesN, icol] .- exDOS[2])])
	end
	push!(temp, temp1)
end

#temp1 tableau longueur unique N, chaque élément se place à un N et prend tous les indices pour chaque éta dans indtempN
