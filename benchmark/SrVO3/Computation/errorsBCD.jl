include("../SrVO3Parameters.jl")

tabN = 3:100
repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64
bz = load_bz(FBZ(), I(d))

prob = DOSProblem(H, float(zero(1.0)), bz)
BCDvalues = zeros(length(Energies), length(tabN))
BCDtimes = zeros(length(Energies), length(tabN))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		temp = @timed AutoBZCore.solve!(cache).value
		BCDvalues[ie, iN] = temp.value
		BCDtimes[ie, iN] = temp.time
	end
	jldsave(repo * "ValuesBCD_N=$(tabN)temp.jld2"; BCDvalues, BCDtimes, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5)
end
jldsave(repo * "ValuesBCD_N=$(tabN).jld2"; BCDvalues, BCDtimes, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5)
