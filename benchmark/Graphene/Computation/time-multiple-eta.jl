include("../GrapheneParameters.jl")

tabN = 3:150
Energies = [0.0, 0.1, 0.5, 0.9, 0.99, 2.0, 2.9, 2.99]
repo = "benchmark/Graphene/Results"

bz = load_bz(FBZ(), I(d))
prob = DOSProblem(H, float(zero(1.0)), bz)
exDOSη = zeros(length(Energies), length(ηlist))
cachetemp = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 401, α = 0.1 / (2π), ΔE = 0.5, η = ηlist[1]))
temp3 = []
@time for (ie, e) in enumerate(Energies)
	cachetemp.domain = e
	temp2 = []
	for (iη, η) in enumerate(ηlist)

		cachetemp.alg = AutoBZCore.BCD(; npt = 401, α = 0.1 / (2π), ΔE = 0.5, η = ηlist[iη])
		exDOSη[ie, iη] = AutoBZCore.solve!(cachetemp).value

	end

	jldsave(repo * "exDOS_Smearing_Multiple.jld2"; exDOSη, α = 0.1 / (2π), ΔE = 0.5, ηlist)
end

#! BCD values
prob = DOSProblem(H, float(zero(1.0)), bz)
BCDtimes = zeros((length(Energies), length(tabN)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = ηlist[1]))
	for (ie, e) in enumerate(Energies)
		#for (iη, η) in enumerate(ηlist)
		cache.alg = AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = ηlist[1])

		cache.domain = e
		BCDtimes[ie, iN] = @elapsed AutoBZCore.solve!(cache).value
		#end
	end
	jldsave(repo * "TimesBCD_N$(tabN).jld2"; BCDtimes = BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = 0, exDOSη)
end
#temperror = [[abs.(BCDvalues[ie, iη, :] .- exDOSη[ie, iη]) for iη in eachindex(ηlist)] for ie in eachindex(Energies)]

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
tolerances = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]))
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
	jldsave(repo * "TestValuesIAI.jld2"; IAIvalues, IAINs, Energies, ηlist, tolerances, exDOSη)
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