include("../SrVO3Parameters.jl")

tabN = 3:300
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62]
repo = "benchmark/SrVO3/Results/"
ηlist = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]), rev = true)
tolerances = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:5]), rev = true)
bz = load_bz(CubicSymIBZ(), I(d))
prob = DOSProblem(H, float(zero(1.0)), bz)


p0 = (; η = 1e-1, ω = 0.01)
greens_function(k, h_k, (; η, ω)) = -imag(tr(inv((ω + im * η) * I - h_k))) / (π * (2π)^d)
prototype = let k = AutoBZCore.FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end

integrand = FourierIntegralFunction(greens_function, H, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-3)

if !isfile(repo * "exDOS_Smearing_Multipleless.jld2")
	exDOSη = zeros(length(Energies), length(ηlist))
	@time for (iη, η) in enumerate(ηlist)
		cachetemp = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 501, α = 0.1 / (2π), ΔE = 0.5, η = η))
		for (ie, e) in enumerate(Energies)
			cachetemp.domain = e
			exDOSη[ie, iη] = AutoBZCore.solve!(cachetemp).value
		end
	end

	jldsave(repo * "exDOS_Smearing_Multiple_less.jld2"; exDOSη, Energies, npt = 501, α = 0.1 / (2π), ΔE = 0.5, ηlist)
else
	exDOSη = load(repo * "exDOS_Smearing_Multiple_less.jld2", "exDOSη")
end
