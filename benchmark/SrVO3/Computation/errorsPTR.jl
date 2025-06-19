include("../SrVO3Parameters.jl")
for (i, arg) in enumerate(ARGS)
	tabN = $arg
end

repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64
bz = load_bz(FBZ(), I(d))
ηlist = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]), rev = true)               # 10 meV (scattering amplitude)

best_of(temp) = temp[argmin([norm(temp[i, :] - exDOS, Inf) for i in axes(temp, 1)]), :]
PTRvalues = zeros(length(Energies), length(ηlist), length(tabN))
PTRtimes = zeros(length(Energies), length(ηlist), length(tabN))
@time for (iN, N) in enumerate(tabN)

	for (iη, η) in enumerate(ηlist)
		for (ie, E) in enumerate(Energies)
			temp = @timed dos_solver_ptr(N, η)(E)
			PTRvalues[ie, iη, iN] = temp.value
			PTRvalues[ie, iη, iN] = temp.time
		end
	end
	jldsave(repo * "PTRValues_N=$(tabN)temp.jld2"; PTRvalues, PTRtimes, Energies, ηlist)
end
jldsave(repo * "PTRValues_N=$(tabN).jld2"; PTRvalues, PTRtimes, Energies, ηlist)