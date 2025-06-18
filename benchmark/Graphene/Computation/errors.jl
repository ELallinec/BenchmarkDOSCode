include("../GrapheneParameters.jl")

tabN = 3:150
Energies = [0.0, 0.1, 0.5, 0.9, 0.99, 2.0, 2.9, 2.99]
exDOS = dos_graphene_exact.(Energies) #Reference Dos
bz = load_bz(FBZ(), I(d))
η = 0
prob = DOSProblem(H, float(zero(1.0)), bz)
BCDvalues = zeros((length(tabN), length(Energies)))

BCDtimes = zeros((length(tabN), length(Energies)))
@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5, η = η))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		temp = @timed AutoBZCore.solve!(cache).value
		BCDvalues[iN, ie] = temp.value
		BCDtimes[iN, ie] = temp.time
	end
	jldsave("benchmark/Graphene/Results/ValuesBCDNoSmearing_N$(tabN)temp.jld2"; BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = η, exDOS)
end
jldsave("benchmark/Graphene/Results/ValuesBCDNoSmearing_N$(tabN)temp.jld2"; BCDvalues, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5, η = η, exDOS)

LTvalues = zeros((length(tabN), length(Energies)))
LTtimes = zeros((length(tabN), length(Energies)))

@time for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.LT(; npt = N))
	for (ie, e) in enumerate(Energies)
		cache.domain = e
		temp = @timed AutoBZCore.solve!(cache).value
		LTvalues[iN, ie] = temp.value
		LTtimes[iN, ie] = temp.time
	end
	jldsave("benchmark/Graphene/Results/ValuesLT_N$(tabN)temp.jld2"; LTvalues, Energies, tabN, exDOS)
end


best_of(temp) = temp[argmin([norm(temp[i, :] - exDOS, Inf) for i in axes(temp, 1)]), :]
PTRvalues = zeros(length(tabN), length(Energies))
PTRtimes = zeros(length(tabN), length(Energies), length(ηlist))
PTRoptimη = zeros(length(tabN), length(Energies))
temp1 = []
@time for n in tabN
	temp = zeros(length(ηlist), length(Energies))

	for (iη, η) in enumerate(ηlist)
		for (ie, E) in enumerate(Energies)
			temp2 = @timed dos_solver_ptr(n, η)(E)
			temp[iη, ie] = temp2.value
			PTRtimes[n-2, ie, iη] = temp2.time
		end
	end
	temp2 = [argmin(abs.(temp[:, ie] .- exDOS[ie])) for ie in eachindex(Energies)]
	PTRoptimη[n-2, :] = ηlist[temp2]

	PTRvalues[n-2, :] = best_of(temp)
	jldsave("benchmark/Graphene/Results/ValuesPTR_N$(tabN)temp.jld2"; PTRvalues, PTRtimes, Energies, tabN, ηlist, exDOS)
end
PTRtime = mean.(eachrow(reduce(hcat, [mean.(eachrow(PTRtimes[:, i, :])) for i in eachindex(Energies)])))
jldsave("benchmark/Graphene/Results/ValuesPTR_N$(tabN)temp.jld2"; PTRvalues, PTRtimes, PTRtime, Energies, tabN, ηlist, exDOS)

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
	jldsave("benchmark/Graphene/Results/ValuesIAItemp.jld2"; IAIvalues, IAINs, IAItimes, Energies, ηlist, tolerances)
end
jldsave("benchmark/Graphene/Results/ValuesIAItemp.jld2"; IAIvalues, IAINs, IAItimes, Energies, ηlist, tolerances)
allerrorsIAI = []
tabNIAI = []
alltabtol = []
alltabη = []
IAItime = []
@time for ie in eachindex(Energies)
	IAIvalue = IAIvalues[ie, :, :]
	IAIN = IAINs[ie, :, :]


	DictIAI = Dict{Float64, Vector{Float64}}()
	Dicttols = Dict{Float64, Vector{Float64}}()
	Dictη = Dict{Float64, Vector{Float64}}()
	Dicttime = Dict{Float64, Vector{Float64}}()
	for (iη, itol) in Iterators.product(eachindex(ηlist), eachindex(tolerances))
		k   = round(IAINs[ie, iη, itol]^(1 / d))
		v   = IAIvalues[ie, iη, itol]
		η   = ηlist[iη]
		tol = tolerances[itol]
		ti  = IAItimes[ie, iη, itol]

		push!(get!(DictIAI, k, Float64[]), v)
		push!(get!(Dicttols, k, Float64[]), tol)
		push!(get!(Dictη, k, Float64[]), η)
		push!(get!(Dicttime, k, Float64[]), ti)
	end
	Dicttols2 = Dict(sort(collect(Dicttols), by = x -> x[1]))
	DictIAI2 = sort(collect(DictIAI), by = x -> x[1])
	Dictη2 = Dict(sort(collect(Dictη), by = x -> x[1]))
	Dicttime2 = Dict(sort(collect(Dicttime), by = x -> x[1]))
	tabNtemp = []
	tabtol = []
	tabη = []
	errorsIAI = []
	IAIti = []
	for (k, v) in DictIAI2
		push!(tabNtemp, k)
		push!(tabη, Dictη2[k][argmin(abs.(v .- exDOS[ie]))])
		push!(tabtol, Dicttols2[k][argmin(abs.(v .- exDOS[ie]))])
		push!(errorsIAI, minimum(abs.(v .- exDOS[ie])))
		push!(IAIti, Dicttime2[k][argmin(abs.(v .- exDOS[ie]))])
	end

	push!(allerrorsIAI, errorsIAI)
	push!(tabNIAI, tabNtemp)
	push!(alltabtol, tabtol)
	push!(alltabη, tabη)
	push!(IAItime, IAIti)
end
jldsave("benchmark/Graphene/Results/ErrorsIAItemp.jld2"; IAIvalues, IAINs, IAItime, Energies, ηlist, tolerances, allerrorsIAI, tabNIAI, alltabtol, alltabη)

