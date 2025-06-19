include("../SrVO3Parameters.jl")
tabN = 3:100
repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64
bz = load_bz(FBZ(), I(d))
ηlist = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]), rev = true)               # 10 meV (scattering amplitude)
tolerances = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:4]), rev = true)
#=

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
jldsave(repo * "ValuesBCD_N=$(tabN).jld2"; BCDvalues, BCDtimes, Energies, tabN)

bz = load_bz(FBZ(), I(d))
prob = DOSProblem(H, float(zero(1.0)), bz)
LTvalues = zeros(length(Energies), length(tabN))
LTtimes = zeros(length(Energies), length(tabN))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.LT(; npt = N))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		temp = @timed AutoBZCore.solve!(cache).value
		LTvalues[ie, iN] = temp.value
		LTtimes[ie, iN] = temp.time
	end
	jldsave(repo * "ValuesLT_N=$(tabN)temp.jld2"; LTvalues, LTtimes, Energies, tabN)
end
jldsave(repo * "ValuesLT_N=$(tabN).jld2"; LTvalues, LTtimes, Energies, tabN)
=#


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

IAIvalues = zeros(SVector{length(Energies), SVector{2, Float64}}, length(ηlist), length(tolerances))
@time for (iη, η) in enumerate(ηlist)
	for (itol, tol) in enumerate(tolerances)
		IAIvalues[iη, itol] = dos_solver_iai(η, tol).(Energies)
	end
end


#=

allN = [getindex.(IAIvalues[i], 2) for i in eachindex(IAIvalues)]
allN2 = [getindex.(allN, i) for i in eachindex(Energies)]
allN23rt = [round.(allN2[i] .^ (1 / 3)) for i in eachindex(allN2)]
allvalues = [getindex.(IAIvalues[i], 1) for i in eachindex(IAIvalues)]
allvalues2 = [getindex.(allvalues, i) for i in eachindex(Energies)]

finalIAIvalues = []
finaltabN = []
fig1, ax = subplots()
for ie in eachindex(Energies)
	temp = Dict{Float64, Vector{Float64}}()
	for (k, v) in zip(allN23rt[ie], allvalues2[ie])
		push!(get!(temp, k, Float64[]), v)
	end
	temp = sort(collect(temp), by = x -> x[1])
	# 2) post‑process so that singletons become Float64, multiples stay Vector

	# → Dict(2.0 => [1.0, 3.0], 5.0 => 2.0, 1.0 => 0.0)
	minima = (
		k => (v isa AbstractVector ? minimum(v) : v)
		for (k, v) in temp
	)

	# 4) Extract keys and min‐values into vectors (in insertion‐order)

	tabN = []
	tabvalues = []
	for (k, v) in temp
		push!(tabN, k)
		push!(tabvalues, minimum(abs.(v .- exDOS[ie])))
	end
	push!(finalIAIvalues, tabvalues)
	push!(finaltabN, tabN)
	ax.plot(tabN, tabvalues, label = "$(Energies[ie])")

end
legend()
yscale(:log)
xscale(:log)


for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(3:50, abs.(BCDvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.scatter(3:50, abs.(LTvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT")
	ax.scatter(3:50, abs.(PTRvalues[:, i] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.scatter(finaltabN[i], abs.(finalIAIvalues[i]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI")
	yscale(:log)
	xlim(3, 50)
	xlabel("N")
	ylabel("Relative error")
	ax.legend()
	savefig("benchmark/SrVO3/fig/RelErrorsVsNSrVO3E=$(Energies[i]).png")
end





# Parameters (tune as needed)
const Nₖ = 50         # grid resolution per dimension
const tol = 1e-2       # energy tolerance for matching eigenvalues
const C = 1.0                # contour-validity constant
#const Energies = [0.01, 0.1, 0.5, 0.9, 0.99, 2.0, 2.9, 2.99]

# Placeholder: user must compute S2 = sum(|R|^2 * ||H_R||) for their model
#const S2 = 3 * a^2 * t  # e.g. S2 = 3 * a^2 * t for SrVO3e with unit a and t

# Initialize storage for minimal gradients
mE = fill(Inf, length(Energies))

# Brillouin zone grid
kx_vals = range(-π, π; length = Nₖ) / 2π
ky_vals = range(-π, π; length = Nₖ) / 2π
kz_vals = range(-π, π; length = Nₖ) / 2π

# Main loop (single-threaded)
@time for kx in kx_vals, ky in ky_vals, kz in kz_vals
	# Get Hamiltonian and its derivatives from user-provided functions
	Hk = H([kx, ky, kz])
	dHkx, dHky, dHkz = dH([kx, ky, kz])[2]

	# Diagonalize H(k)
	evals, evecs = eigen(Hk)

	# For each target energy, find close eigenvalues and update gradient
	for (i, E0) in enumerate(Energies)
		for (j, λ) in enumerate(evals)
			if abs(real(λ) - E0) < tol
				u = evecs[:, j]
				# Hellmann–Feynman gradient components
				dλ_dkx = real(dot(u, dHkx * u))
				dλ_dky = real(dot(u, dHky * u))
				dλ_dkz = real(dot(u, dHkz * u))
				grad = sqrt(dλ_dkx^2 + dλ_dky^2 + dλ_dkz^2)
				# Update minimum gradient
				if grad < mE[i]
					mE[i] = grad
				end
			end
		end
	end
end

# Compute and print optimal ΔE* and α*
println("Energy, m_E, ΔE*, α*")
for (E0, m) in zip(Energies, mE)
	if isfinite(m) && m > 0
		ΔE = C + S2 / m
		α = ΔE / S2
		println("E0=", E0, " m=", m, " ΔE=", ΔE, " α=", α)
	else
		println("%6.2f, no match or zero gradient\n", E0)
	end
end
=#