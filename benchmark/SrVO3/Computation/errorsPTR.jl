include("../SrVO3Parameters.jl")

tabN = 201:400
repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64

ηlist = sort(logrange(1e-4, 1, 14), rev = true)
bz = load_bz(CubicSymIBZ(), "data/svo.wout")
p0 = (; η = 1e-2, ω = 12.44) # initial parameters
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end
integrand = FourierIntegralFunction(greens_function, H, prototype)

function dos_solver_ptr( η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0), AutoPTR(; a = 1e-5))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end
best_of(temp) = temp[argmin([norm(temp[i, :] - exDOS, Inf) for i in axes(temp, 1)]), :]
PTRvalues = zeros(length(Energies), length(ηlist), length(tabN))
PTRtimes = zeros(length(Energies), length(ηlist), length(tabN))
for (iN, N) in enumerate(tabN)
	@time for (ie, E) in enumerate(Energies)
		for (iη, η) in enumerate(ηlist)

			temp = @timed dos_solver_ptr(N, η)(E)
			PTRvalues[ie, iη, iN] = -imag(temp.value) / det(bz.B) / π
			PTRtimes[ie, iη, iN] = temp.time
		end
	end
	jldsave(repo * "PTRValues_N=$(tabN)_final.jld2"; PTRvalues, PTRtimes, Energies, ηlist)
end
#jldsave(repo * "PTRValues_N=$(tabN).jld2"; PTRvalues, PTRtimes, Energies, ηlist)