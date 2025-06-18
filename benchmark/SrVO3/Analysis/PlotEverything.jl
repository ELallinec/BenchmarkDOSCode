using JLD2, LinearAlgebra, PyPlot, AutoBZCore
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 18)
PyPlot.rc("figure", figsize = (9, 6 * 9 / 8))
tabindex = vcat(1:77, 79:101)

@load "benchmark/SrVO3/Results/RefDos.jld2"
@load "benchmark/SrVO3/Results/BCDerrors+time.jld2"
@load "benchmark/SrVO3/Results/LTerrors+time.jld2"
@load "benchmark/SrVO3/Results/PTRerrors+time.jld2"
@load "benchmark/SrVO3/Results/IAI.jld2"


temp = Matrix.(reduce.(hcat, IAIvalues))
for i in CartesianIndices(temp)

	temp[i][2, :] = Int64.(ceil.(temp[i][2, :] .^ (1 / 3)))

end

tabmaxNfull = zeros(size(temp))
tabmeanNfull = zeros(size(temp))
taberrorsfull = zeros(size(temp))


for i in CartesianIndices(temp)

	tabmaxNfull[i] = Int64(maximum(temp[i][2, :]))
	tabmeanNfull[i] = Int64(ceil(mean(temp[i][2, :])))
	taberrorsfull[i] = norm(temp[i][1, :] - RefDos, Inf)

end
tabNmaxfull = sort(unique(tabmaxNfull))
tabNmeanfull = sort(unique(tabmeanNfull))
allindicesmaxNfull = [findall(i -> i == unique(tabNmaxfull)[j], tabmaxNfull) for j in CartesianIndices(tabNmaxfull)]
allindicesmeanNfull = [findall(i -> i == unique(tabNmeanfull)[j], tabmeanNfull) for j in CartesianIndices(tabNmeanfull)]
errorsmaxNIAIfull = zeros(size(tabNmaxfull))
errorsmeanNIAIfull = zeros(size(tabNmeanfull))
optimηmaxNIAIfull = zeros(size(tabNmaxfull))
optimηmeanNIAIfull = zeros(size(tabNmeanfull))
optimtolmaxNIAIfull = zeros(size(tabNmaxfull))
optimtolmeanNIAIfull = zeros(size(tabNmeanfull))
for (iv, vect) in enumerate(allindicesmaxNfull)
	errorsmaxNIAIfull[iv] = minimum(taberrorsfull[vect])
	#println(argmin(taberrorsfull[vect]))
	optimηmaxNIAIfull[iv] = ηlist[vect[argmin(taberrorsfull[vect])][2]]
	optimtolmaxNIAIfull[iv] = tolerances[vect[argmin(taberrorsfull[vect])][1]]
end

for (iv, vect) in enumerate(allindicesmeanNfull)
	errorsmeanNIAIfull[iv] = minimum(taberrorsfull[vect])
	optimηmeanNIAIfull[iv] = ηlist[vect[argmin(taberrorsfull[vect])][2]]
	optimtolmeanNIAIfull[iv] = tolerances[vect[argmin(taberrorsfull[vect])][1]]
end


###


tabmaxNnosing = zeros(size(temp))
tabmeanNnosing = zeros(size(temp))
taberrorsnosing = zeros(size(temp))


for i in CartesianIndices(temp)

	tabmaxNnosing[i] = Int64(maximum(temp[i][2, tabindex]))
	tabmeanNnosing[i] = Int64(ceil(mean(temp[i][2, tabindex])))
	taberrorsnosing[i] = norm(temp[i][1, tabindex] - RefDos[tabindex], Inf)

end
tabNmaxnosing = sort(unique(tabmaxNnosing))
tabNmeannosing = sort(unique(tabmeanNnosing))
allindicesmaxNnosing = [findall(i -> i == unique(tabNmaxnosing)[j], tabmaxNnosing) for j in CartesianIndices(tabNmaxnosing)]
allindicesmeanNnosing = [findall(i -> i == unique(tabNmeannosing)[j], tabmeanNnosing) for j in CartesianIndices(tabNmeannosing)]
errorsmaxNIAInosing = zeros(size(tabNmaxnosing))
errorsmeanNIAInosing = zeros(size(tabNmeannosing))
optimηmaxNIAInosing = zeros(size(tabNmaxnosing))
optimηmeanNIAInosing = zeros(size(tabNmeannosing))
optimtolmaxNIAInosing = zeros(size(tabNmaxnosing))
optimtolmeanNIAInosing = zeros(size(tabNmeannosing))
for (iv, vect) in enumerate(allindicesmaxNnosing)
	errorsmaxNIAInosing[iv] = minimum(taberrorsnosing[vect])
	optimηmaxNIAInosing[iv] = ηlist[vect[argmin(taberrorsnosing[vect])][2]]
	optimtolmaxNIAInosing[iv] = ηlist[vect[argmin(taberrorsnosing[vect])][1]]
end

for (iv, vect) in enumerate(allindicesmeanNnosing)
	errorsmeanNIAInosing[iv] = minimum(taberrorsnosing[vect])
	optimηmeanNIAInosing[iv] = ηlist[vect[argmin(taberrorsnosing[vect])][2]]
	optimtolmeanNIAInosing[iv] = ηlist[vect[argmin(taberrorsnosing[vect])][1]]
end
