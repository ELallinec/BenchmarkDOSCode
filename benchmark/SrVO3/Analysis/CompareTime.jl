include("../SrVO3Parameters.jl")

using JLD2
path = "benchmark/SrVO3/ResultsNonSmearing/"
d = 3
exDOS, Energies = load(path * "exDOS_N=501.jld2", "exDOS", "Energies")
BCDvalues, BCDtimes = load(path * "ValuesBCD_N=3:200_final.jld2", "BCDvalues", "BCDtimes")
LTvalues, LTtimes = load(path * "ValuesLT_N=3:200_final.jld2", "LTvalues", "LTtimes")
IAIvalues, IAItimes, IAINs, tolerances, ηlist = load(path * "ValuesIAI_final2.jld2", "IAIvalues", "IAItimes", "IAINs", "tolerances", "ηlist")
PTRvalues, PTRtimes = load(path * "ValuesPTR_N=3:400_final.jld2", "PTRvalues", "PTRtimes")

timePTR = [minimum(PTRtimes[:, :, n]) for n in 1:198]
timeBCD = [minimum(BCDtimes[:, n]) for n in 1:198]
timeLT = [minimum(LTtimes[:, n]) for n in 1:198]
PTR_IAI = []
BCD_IAI = []
LT_IAI = []
for i in 1:6
	for (iN, N) in enumerate(tabNIAI[i])
		if N <= 200
			push!(PTR_IAI, timePTR[Int64(N - 2)] / IAItime[i][iN])
			push!(LT_IAI, timeLT[Int64(N - 2)] / IAItime[i][iN])
			push!(BCD_IAI, timeBCD[Int64(N - 2)] / IAItime[i][iN])
		end
	end
end

