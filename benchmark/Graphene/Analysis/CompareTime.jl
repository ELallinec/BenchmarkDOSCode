include("../GrapheneParameters.jl")

@load "benchmark/Graphene/ResultsNonSmearing/ErrorsIAI_final.jld2"
#@load "benchmark/Graphene/Results/ErrorsIAI.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesLT_N3:400_final.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesBCD_N3:400_final.jld2"
@load "benchmark/Graphene/ResultsNonSmearing/ValuesPTR_N3:400_final.jld2"


timePTR = [minimum(PTRtimes[:, :, n]) for n in 1:398]
timeBCD = [minimum(BCDtimes[:, n]) for n in 1:398]
timeLT = [minimum(LTtimes[:, n]) for n in 1:398]
PTR_IAI = []
BCD_IAI = []
LT_IAI = []
for i in 1:8
	for (iN, N) in enumerate(tabNIAI[i])
		if N <= 400
			push!(PTR_IAI, timePTR[Int64(N - 2)] / IAItime[i][iN])
			push!(LT_IAI, timeLT[Int64(N - 2)] / IAItime[i][iN])
			push!(BCD_IAI, timeBCD[Int64(N - 2)] / IAItime[i][iN])
		end
	end
end

