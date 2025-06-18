include("../Mono1DParameters.jl")

#! Tunable parameters
tabN = 3:1500
Energies = [0.0, 1.0, 1.5, 1.9, 1.99]
repo = "benchmark/Mono1D/Results/"

#! Reference DOS
exDOS = @. 1 / (π * sqrt(4 - Energies^2)); #Reference Dos

#! Computations

# BCD values
prob = DOSProblem(H, float(zero(1.0)), bz)
BCDvalues = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = 0))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		BCDvalues[iN, ie] = AutoBZCore.solve!(cache).value
	end
	jldsave(repo * "ValuesBCDNoSmearing_N$(tabN).jld2"; BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOS)
end

#! LT values
LTvalues = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.LT(; npt = N))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		LTvalues[iN, ie] = AutoBZCore.solve!(cache).value
	end
	jldsave(repo * "ValuesLTNoSmearing_N$(tabN).jld2"; LTvalues, Energies, tabN, exDOS)
end

#! PTR values
best_of(temp) = temp[argmin([norm(temp[i, :] - exDOS, Inf) for i in axes(temp, 1)]), :]
PTRvalues = zeros(length(tabN), length(Energies))
PTRoptimη = zeros(length(tabN), length(Energies))
temp1 = []
@time for n in tabN
	temp = zeros(length(ηlist), length(Energies))

	for (iη, η) in enumerate(ηlist)
		for (ie, E) in enumerate(Energies)
			temp[iη, ie] = dos_solver_ptr(n, η)(E)
		end
	end
	temp2 = [argmin(abs.(temp[:, ie] .- exDOS[ie])) for ie in eachindex(Energies)]
	#println(temp2)
	PTRoptimη[n-2, :] = ηlist[temp2]

	PTRvalues[n-2, :] = best_of(temp)
	jldsave(repo * "ValuesPTR_N$(tabN).jld2"; PTRvalues, Energies, tabN, ηlist, exDOS)
end

#! IAI values
IAIvalues = zeros(length(Energies), length(ηlist), length(tolerances))
IAINs = zeros(length(Energies), length(ηlist), length(tolerances))
@time for (itol, tol) in enumerate(tolerances)
	for (iη, η) in enumerate(ηlist)
		for (iE, E) in enumerate(Energies)
			temp = dos_solver_iai(η, tol)(E)
			IAIvalues[iE, iη, itol] = temp[1]
			IAINs[iE, iη, itol] = temp[2]
		end
	end
	jldsave(repo * "ValuesIAI.jld2"; IAIvalues, IAINs, Energies, ηlist, tolerances, exDOS)
end

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
jldsave(repo * "ErrorsIAI.jld2"; IAIvalues, IAINs, Energies, ηlist, tolerances, allerrorsIAI, tabNIAI, alltabtol, alltabη, exDOS)

