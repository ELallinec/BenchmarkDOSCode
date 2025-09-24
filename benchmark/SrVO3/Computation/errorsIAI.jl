include("../SrVO3Parameters.jl")

repo = "benchmark/SrVO3/Results/"
Energies = [11.55, 12.44, 13.19, 13.29, 13.46, 13.62] #Singularities at 11.57, 13.31, 13.64

ηlist = sort(logrange(10^-4, 1, 14), rev = true)
tolerances = sort(unique([i * 10.0^(-j) for i in [1, 5, 10] for j in 1:5]), rev = true)
εF = 12.3958
bz = load_bz(CubicSymIBZ(), "data/svo.wout")
p0 = (; η = 1e-1, ω = εF) # initial parameters
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end
integrand = FourierIntegralFunction(greens_function, H, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-5)


function dos_solver_iai(η, atol)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz1, p0; abstol = atol), EvalCounter(IAI(AuxQuadGKJL(order = 4))))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)
		@SVector [getproperty(temp, :value), getproperty(getproperty(temp, :stats), :numevals)]
	end
end

function dos_solver_ptr(η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0), AutoPTR(; nmin = 1); abstol = 1e-5)
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end

function dos_solver_ptr(η, N)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz1, p0), PTR(npt = N))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end

#=
IAIvalues = zeros(length(Energies), length(ηlist))
IAItimes = zeros(length(Energies), length(ηlist))
IAINs = zeros(length(Energies), length(ηlist))
=#
IAIval = zeros(length(Energies), length(ηlist))
IAItim = zeros(length(Energies), length(ηlist))
IAIN = zeros(length(Energies), length(ηlist))
for (ie, E) in enumerate(Energies)
	@time for (iη, η) in enumerate(ηlist)
		for (itol, tol) in enumerate(tolerances[end:end])
			temp = @timed dos_solver_iai(η, tol)(E)
			IAIval[ie, iη] = -imag(temp.value[1]) / det(bz.B) / π
			IAIN[ie, iη] = temp.value[2]
			IAItim[ie, iη] = temp.time
			jldsave("benchmark/SrVO3/Results/ValuesIAI_Multipleη_1e-5.jld2"; PTRval, PTRtim, ηlist, εF)
		end
	end
end
jldsave("benchmark/SrVO3/Results/ValuesIAI_Multipleη_1e-5.jld2"; IAIval, IAItim, IAIN, ηlist, εF)

PTRval = zeros(length(Energies), length(ηlist))
PTRtim = zeros(length(Energies), length(ηlist))
for (iη, η) in enumerate(ηlist)
	@time for (ie, E) in enumerate(Energies)

		temp = @timed dos_solver_ptr(η)(E)
		PTRval[ie, iη] = -imag(temp.value) / det(bz.B) / π
		PTRtim[ie, iη] = temp.time
		jldsave("benchmark/SrVO3/Results/ValuesPTR_Multipleη_1e-5.jld2"; PTRval, PTRtim, ηlist, εF)
	end
end
#end
#jldsave("benchmark/SrVO3/Results/ValuesIAI_Multipleη_1e-5.jld2"; IAIvalues, IAINs, IAItimes, ηlist, tolerances, Energies)


plot(ηlist, IAItim, marker = :^, label = "IAI", lw = 5, markersize = 20)
plot(ηlist, PTRtim, marker = :v, label = "PTR", lw = 5, markersize = 20)
plot(ηlist[3:7], (ηlist[3:7] .^ -3) / 1000, "--k", label = "lin", lw = 5)
plot(ηlist[8:end], log.(ηlist[8:end] .^ -1) .^ 3 / 10, "--k", label = "log3", lw = 5)
legend()
xlabel("η (eV)")
ylabel("Time (s)")
xscale(:log)
yscale(:log)
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