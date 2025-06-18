include("../Mono1DParameters.jl")

# Tunable parameters
η = 1e-3
tabN = 3:3165
Energies = [0.0, 1.0, 1.5, 1.9, 1.99]
repo = "benchmark/Mono1D/Results/"

# Reference DOS
if !isfile(repo * "exDOS_Smearing.jld2")

	cachetemp = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 1e8 + 1, α = 0.1 / (2π), ΔE = 0.5, η = η))
	exDOSη = []
	@time for (ie, e) in enumerate(Energies)
		cachetemp.domain = e
		push!(exDOSη, AutoBZCore.solve!(cachetemp).value)
	end

	jldsave(repo * "exDOS_Smearing.jld2"; exDOSη, Energies, npt = 1e8 + 1, α = 0.1 / (2π), ΔE = 0.5, η = η)
else
	exDOSη = load(repo * "exDOS_Smearing.jld2", "exDOSη")
end
# Computations


#? BCD values
BCDvalues = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = η))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		BCDvalues[iN, ie] = AutoBZCore.solve!(cache).value
	end
	jldsave(repo * "ValuesBCD_Smearing_N$(tabN).jld2"; BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = η, exDOSη)

end


#? PTR values
PTRvalues = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	for (ie, E) in enumerate(Energies)
		PTRvalues[iN, ie] = dos_solver_ptr(N, η)(E)
	end
	jldsave(repo * "ValuesPTR_Smearing_N$(tabN).jld2"; PTRvalues, Energies, tabN, η = η, exDOSη)
end


#? IAI values
newtolerances = sort(unique([i * 10.0^(-j) for i in 1:0.01:10 for j in 0:14]), rev = true)
#newtolerances = tolerances
IAIvalues = zeros(length(Energies), length(newtolerances))
IAINs = zeros(length(Energies), length(newtolerances))
@time for (itol, tol) in enumerate(newtolerances)
	for (iE, E) in enumerate(Energies)
		result = dos_solver_iai(η, tol)(E)
		IAIvalues[iE, itol] = result[1]
		IAINs[iE, itol] = result[2]
	end
	jldsave(repo * "ValuesIAI_Smearing.jld2"; IAIvalues, IAINs, Energies, η = η, tolerances = newtolerances, exDOSη)
end

#? IAI absolute errors
allerrorsIAI = []
tabNIAI = []
alltabtol = []
for ie in eachindex(Energies)
	IAIvalue = IAIvalues[ie, :]
	IAIN = IAINs[ie, :]


	DictIAI = Dict{Float64, Vector{Float64}}()
	Dicttols = Dict{Float64, Vector{Float64}}()
	for (k, v) in zip(IAIN, IAIvalue)
		push!(get!(DictIAI, k, Float64[]), v)
	end
	for (k, v) in zip(IAIN, newtolerances)
		push!(get!(Dicttols, k, Float64[]), v)
	end
	Dicttols2 = sort(collect(Dicttols), by = x -> x[1])
	DictIAI2 = sort(collect(DictIAI), by = x -> x[1])

	tabNtemp = []
	tabtol = []
	errorsIAI = []
	for (k, v) in DictIAI2
		push!(tabNtemp, k)
		push!(tabtol, Dicttols[k][argmin(abs.(v .- exDOSη[ie]))])
		push!(errorsIAI, minimum(abs.(v .- exDOSη[ie])))
	end

	push!(allerrorsIAI, errorsIAI)
	push!(tabNIAI, tabNtemp)
	push!(alltabtol, tabtol)
end

jldsave(repo * "ErrorsIAI_Smearing.jld2"; IAIvalues, IAINs, Energies, η = η, tolerances = newtolerances, allerrorsIAI, tabNIAI, alltabtol, exDOSη)
