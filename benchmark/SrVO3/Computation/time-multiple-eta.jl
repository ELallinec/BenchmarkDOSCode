include("../SrVO3Parameters.jl")

tabN = 3:100
Energies = [11.55, 12.44, 13.29, 13.62]
repo = "benchmark/SrVO3/Results/"

bz = load_bz(CubicSymIBZ(), I(d))
p0 = (; η = 1e-1, ω = 0.01)
greens_function(k, h_k, (; η, ω)) = -imag(tr(inv((ω + im * η) * I - h_k))) / (π * (2π)^d)
prototype = let k = AutoBZCore.FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end

integrand = FourierIntegralFunction(greens_function, H, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-3)

function dos_solver_iai(η, atol)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = atol), EvalCounter(IAI()))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		@SVector [getproperty(temp, :value), getproperty(getproperty(temp, :stats), :numevals)]
	end
end

function dos_solver_ptr(N, η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0), PTR(npt = N))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end

prob = DOSProblem(H, float(zero(1.0)), bz)
exDOSη = zeros(length(Energies), length(ηlist))
cachetemp = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 201, α = 0.1 / (2π), ΔE = 0.5, η = 0))
temp3 = []
@time for (ie, e) in enumerate(Energies)
	cachetemp.domain = e
	temp2 = []
	for (iη, η) in enumerate(ηlist)

		cachetemp.alg = AutoBZCore.BCD(; npt = 401, α = 0.1 / (2π), ΔE = 0.5, η = ηlist[iη])
		exDOSη[ie, iη] = AutoBZCore.solve!(cachetemp).value

	end

	jldsave(repo * "exDOS_Smearing_Multipletemp.jld2"; exDOSη, α = 0.1 / (2π), ΔE = 0.5, ηlist)
end
jldsave(repo * "exDOS_Smearing_Multiple.jld2"; exDOSη, α = 0.1 / (2π), ΔE = 0.5, ηlist)

exDOS = zeros(length(Energies))
@time for (ie, E) in enumerate(Energies)
	cachetemp.domain = E
	cachetemp.alg = AutoBZCore.BCD(; npt = 401, α = 0.1 / (2π), ΔE = 0.5, η = 0)
	exDOS[ie] = AutoBZCore.solve!(cachetemp).value
end
#! BCD values
prob = DOSProblem(H, float(zero(1.0)), bz)
BCDvalues = zeros((length(Energies), length(ηlist), length(tabN)))
BCDtimes = zeros((length(Energies), length(ηlist), length(tabN)))
@time for (iN, N) in enumerate(tabN)
	for (iη, η) in enumerate(ηlist)
		cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = η))
		for (ie, e) in enumerate(Energies)

			cache.alg = AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = η)

			cache.domain = e
			temp = @timed AutoBZCore.solve!(cache).value
			BCDvalues[ie, iη, iN] = temp.value
			BCDtimes[ie, iη, iN] = temp.time

		end
	end
	jldsave(repo * "ValuesBCD_Multipleη_N$(tabN)temp.jld2"; BCDtimes = BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOSη)
end
jldsave(repo * "ValuesBCD_Multipleη_N$(tabN).jld2"; BCDtimes, BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOSη)
#temperror = [[abs.(BCDvalues[ie, iη, :] .- exDOSη[ie, iη]) for iη in eachindex(ηlist)] for ie in eachindex(Energies)]

#! PTR values
best_of(temp) = temp[argmin([norm(temp[i, :] - exDOS, Inf) for i in axes(temp, 1)]), :]
PTRvalues = zeros(length(Energies), length(ηlist), length(tabN))
PTRtimes = zeros(length(Energies), length(ηlist), length(tabN))

@time for n in tabN
	for (iη, η) in enumerate(ηlist)
		for (ie, E) in enumerate(Energies)
			ttemp = @timed dos_solver_ptr(n, η)(E)
			PTRvalues[ie, iη, n-2] = ttemp.value
			PTRtimes[ie, iη, n-2] = ttemp.time
		end
	end
	jldsave(repo * "ValuesPTR_Multipleη_N$(tabN)temp.jld2"; PTRtimes, PTRvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOSη)
end
jldsave(repo * "ValuesPTR_Multipleη_N$(tabN).jld2"; PTRtimes, PTRvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOSη)
#! IAI values
IAIvalues = zeros(length(Energies), length(ηlist), length(tolerances))
IAINs = zeros(length(Energies), length(ηlist), length(tolerances))
IAItimes = zeros(length(Energies), length(ηlist), length(tolerances))
@time for (itol, tol) in enumerate(tolerances)
	for (iη, η) in enumerate(ηlist)
		for (iE, E) in enumerate(Energies)
			temp = @timed dos_solver_iai(η, tol)(E)
			IAIvalues[iE, iη, itol] = temp.value[1]
			IAINs[iE, iη, itol] = temp.value[2]
			IAItimes[iE, iη, itol] = temp.time
		end
	end
	jldsave(repo * "ValuesIAI_Multipleηtemp.jld2"; IAIvalues, IAINs, IAItimes, Energies, ηlist, tolerances, exDOSη)
end
jldsave(repo * "ValuesIAI_Multipleη.jld2"; IAIvalues, IAINs, IAItimes, Energies, ηlist, tolerances, exDOSη)

#! IAI absolute error
allerrorsIAI = []
tabNIAI = []
alltabtol = []
alltabη = []
@time for ie in eachindex(Energies)
	IAIvalue = IAIvalues[ie, :, :]
	IAIN = IAINs[ie, :, :]


	DictIAI = Dict{Float64, Vector{Float64}}()
	Dicttols = Dict{Float64, Vector{Float64}}()
	Dictη = Dict{Float64, Vector{Float64}}()

	for (iη, itol) in Iterators.product(eachindex(ηlist), eachindex(tolerances))
		k   = round(IAINs[ie, iη, itol]^(1 / d))
		v   = IAIvalues[ie, iη, itol]
		η   = ηlist[iη]
		tol = tolerances[itol]

		push!(get!(DictIAI, k, Float64[]), v)
		push!(get!(Dicttols, k, Float64[]), tol)
		push!(get!(Dictη, k, Float64[]), η)
	end
	Dicttols2 = Dict(sort(collect(Dicttols), by = x -> x[1]))
	DictIAI2 = sort(collect(DictIAI), by = x -> x[1])
	Dictη2 = Dict(sort(collect(Dictη), by = x -> x[1]))
	tabNtemp = []
	tabtol = []
	tabη = []
	errorsIAI = []
	for (k, v) in DictIAI2
		push!(tabNtemp, k)
		push!(tabη, Dictη2[k][argmin(abs.(v .- exDOS[ie]))])
		push!(tabtol, Dicttols2[k][argmin(abs.(v .- exDOS[ie]))])
		push!(errorsIAI, minimum(abs.(v .- exDOS[ie])))
	end

	push!(allerrorsIAI, errorsIAI)
	push!(tabNIAI, tabNtemp)
	push!(alltabtol, tabtol)
	push!(alltabη, tabη)
end
@load "benchmark/SrVO3/Results/ValuesBCD_Multipleη_N3:100.jld2"
@load "benchmark/SrVO3/Results/ValuesIAI_Multipleη.jld2"
@load "benchmark/SrVO3/Results/ValuesPTR_Multipleη_N3:100.jld2"
function mini(list)
	if isempty(list)
		return Inf
	else
		return minimum(list)
	end
end

threshold = 1e-2
indPTR = []

for ie in eachindex(Energies)
	temp1 = []
	for iη in eachindex(ηlist)
		temp2 = findall(ind -> abs(PTRvalues[ie, iη, ind] - exDOS[ie]) / abs(exDOS[ie]) <= threshold, eachindex(PTRvalues[ie, iη, :]))
		push!(temp1, temp2)
	end
	push!(indPTR, temp1)
end


tPTR = [[(PTRtimes[i, j, indPTR[i][j]]) for j in eachindex(ηlist)] for i in eachindex(Energies)]
mintimesPTR = [[mini(tPTR[i][j]) for j in eachindex(ηlist)] for i in eachindex(Energies)]


indIAI = []
for ie in eachindex(Energies)
	temp1 = []
	for iη in eachindex(ηlist)
		temp2 = findall(ind -> abs(IAIvalues[ie, iη, ind] - exDOS[ie]) / abs(exDOS[ie]) <= threshold, eachindex(IAIvalues[ie, iη, :]))
		push!(temp1, temp2)
	end
	push!(indIAI, temp1)
end


tIAI = [[(IAItimes[i, j, indIAI[i][j]]) for j in eachindex(ηlist)] for i in eachindex(Energies)]
mintimesIAI = [[mini(tIAI[i][j]) for j in eachindex(ηlist)] for i in eachindex(Energies)]

indBCD = []
for ie in eachindex(Energies)
	temp1 = []
	for iη in eachindex(ηlist)
		temp2 = findall(ind -> abs(BCDvalues[ie, iη, ind] - exDOS[ie]) / abs(exDOS[ie]) <= threshold, eachindex(BCDvalues[ie, iη, :]))
		push!(temp1, temp2)
	end
	push!(indBCD, temp1)
end


tBCD = [[(BCDtimes[i, j, indBCD[i][j]]) for j in eachindex(ηlist)] for i in eachindex(Energies)]
mintimesBCD = [[mini(tBCD[i][j]) for j in eachindex(ηlist)] for i in eachindex(Energies)]


for (ie, E) in enumerate(Energies)
	fig, ax = subplots()
	ax.scatter(ηlist, mintimesBCD[ie], marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD")
	ax.plot(ηlist, mintimesBCD[ie], marker = "v", color = (204 / 255, 121 / 255, 167 / 255))
	ax.scatter(ηlist, mintimesPTR[ie], marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR")
	ax.plot(ηlist, mintimesPTR[ie], marker = "+", color = (0, 158 / 255, 115 / 255))
	ax.scatter(ηlist, mintimesIAI[ie], marker = "^", color = :black, label = "IAI")
	ax.plot(ηlist, mintimesIAI[ie], marker = "^", color = :black)
	ax.set_xscale(:log)
	ax.set_yscale(:log)
	ax.legend()
	ax.set_xlabel("η (eV)")
	ax.set_ylabel("Wall-Clock time (s)")
	ax.set_title("E=$(E)")
	#savefig(repo * "../fig/Wall-ClockTimevsEtaE=$(E)tol=$(threshold).png")
end
