include("../SrVO3Parameters.jl")

tabN = 3:100
repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64
bz = load_bz(CubicSymIBZ(), I(d))

prob = DOSProblem(H, float(zero(1.0)), bz)
BCDvalues = zeros(length(Energies), length(tabN))
BCDtimes = zeros(length(Energies), length(tabN))

BCDval = zeros(length(Energies), length(ηlist))
BCDtim = zeros(length(Energies), length(ηlist))
for (ie, e) in enumerate(Energies)

	@time for (iη, η) in enumerate(ηlist)
		cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 100', α = 0.1 / (2π), ΔE = 0.5, η = η))
		cache.domain = e
		temp1 = @btimed AutoBZCore.solve!($cache).value
		BCDval[ie, iη] = temp1.value
		BCDtim[ie, iη] = temp1.time
	end
	#jldsave(repo * "ValuesBCD_Multipleη_1e-5.jld2"; BCDval, BCDtim, Energies, α = 0.1 / (2π), ΔE = 0.5, ηlist)
end
#jldsave(repo * "ValuesBCD_N=$(tabN).jld2"; BCDvalues, BCDtimes, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5)
