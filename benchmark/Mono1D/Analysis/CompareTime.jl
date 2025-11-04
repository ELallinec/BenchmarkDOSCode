include("../Mono1DParameters.jl")


@load "benchmark/Mono1D/Results/ValuesBCD_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesPTR_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesLT_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesIAI_final.jld2"

timePTR = [minimum(PTRtimes[:, :, n]) for n in 1:1498]
timeBCD = [minimum(BCDtimes[:, n]) for n in 1:1498]
timeLT = [minimum(LTtimes[:, n]) for n in 1:1498]
PTR_IAI = []
BCD_IAI = []
LT_IAI = []
for i in 1:5
	for (iN, N) in enumerate(tabNIAI[i])
		if N <= 1500
			push!(PTR_IAI, timePTR[Int64(N - 2)] / IAItime[i][iN])
			push!(LT_IAI, timeLT[Int64(N - 2)] / IAItime[i][iN])
			push!(BCD_IAI, timeBCD[Int64(N - 2)] / IAItime[i][iN])
		end
	end
end

