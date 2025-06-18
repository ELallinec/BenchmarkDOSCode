include("../Mono1DParameters.jl")

# Tunable Parameters
tabN = 3:1500
Energies = [0.0, 1.0, 1.5, 1.9, 1.99]
repo = "benchmark/Mono1D/Results/"
# Reference DOS
exDOS = @. 1 / (π * sqrt(4 - Energies^2));

# COmputation 


#? BCD times
prob = DOSProblem(H, float(zero(1.0)), bz)
BCDtimes = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = 0))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		BCDtimes[iN, ie] = @elapsed AutoBZCore.solve!(cache).value
	end
	jldsave(repo * "timesBCD_N$(tabN).jld2"; BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOS)
end


#? LT times
LTtimes = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.LT(; npt = N))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		LTtimes[iN, ie] = @elapsed AutoBZCore.solve!(cache).value
	end
	jldsave(repo * "timesLT_N$(tabN).jld2"; LTvalues, Energies, tabN, exDOS)
end


#? PTR times
PTRtimes = zeros(length(tabN), length(Energies))
@time for n in tabN
	temp = zeros(length(ηlist), length(Energies))

	for (ie, E) in enumerate(Energies)
		PTRtimes[n-2, ie] = @elapsed dos_solver_ptr(n, η)(E)
	end
	jldsave(repo * "TimesPTR_N$(tabN).jld2"; PTRvalues, Energies, tabN, ηlist, exDOS)
end


#? IAI times (need precomputation of optimal η and tolerance)
@load repo * "ErrorsIAI.jld2"
IAItimes = [zeros(length(tabNIAI)) for ie in eachindex(Energies)]

for (ie, E) in enumerate(Energies)
	for (iN, N) in enumerate(tabNIAI)
		IAItimes[ie][iN] = @elapsed dos_solver_iai(alltabtol[ie, iN], alltabη[ie, iN])
	end
	jldsave(repo * "TimesIAI.jld2"; IAItimes, Energies, alltabtol, talltabη, exDOS)
end
