include("../SrVO3Parameters.jl")

Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62]
repo = "benchmark/SrVO3/Results/"
bz = load_bz(CubicSymIBZ(), I(d))
prob = DOSProblem(H, float(zero(1.0)), bz)
npt = 1001
if !isfile(repo * "exDOS_Smearing_Multiple.jld2")
	exDOS = zeros(length(Energies))
	cachetemp = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = npt, α = 0.1 / (2π), ΔE = 0.5))
	@time for (ie, e) in enumerate(Energies)
		cachetemp.domain = e
		exDOS[ie] = AutoBZCore.solve!(cachetemp).value
	end

	jldsave(repo * "exDOS_N=$(npt).jld2"; exDOS, Energies, npt = npt, α = 0.1 / (2π), ΔE = 0.5)
else
	exDOS = load(repo * "exDOS_N=$(npt).jld2", "exDOS")
end

Energies = [11.569]
bz1 = load_bz(CubicSymIBZ(), I(d))
prob1 = DOSProblem(H, float(zero(1.0)), bz1)
@time cachetemp = AutoBZCore.init(prob1, AutoBZCore.BCD(; npt = 701, α = 0.1 / (2π), ΔE = 0.5, η = 0));
exDOS = zeros(length(Energies))
@time for (ie, e) in enumerate(Energies)
	cachetemp.domain = e
	exDOS[ie] = AutoBZCore.solve!(cachetemp).value
end