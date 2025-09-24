include("../GrapheneParameters.jl")

tabN = 3:400
Energies = [0.0, 0.1, 0.5, 0.9, 0.99, 2.0, 2.9, 2.99]
exDOS = dos_graphene_exact.(Energies) #Reference Dos
bz = load_bz(FBZ(), I(d))
prob = DOSProblem(H, 0.99, bz)

BCDvalues = zeros(length(Energies), length(tabN))
BCDtimes = zeros(length(Energies), length(tabN))
for (iN, N) in enumerate(tabN)
	cache = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = N, α = 0.1 / (2π), ΔE = 0.5))
	@time for (ie, e) in enumerate(Energies)
		cache.domain = e
		temp = @timed AutoBZCore.solve!(cache).value
		BCDvalues[ie, iN] = temp.value
		BCDtimes[ie, iN] = temp.time
	end
	jldsave("benchmark/Graphene/ResultsNonSmearing/ValuesBCD_N$(tabN)_final2.jld2"; BCDvalues, BCDtimes, Energies, tabN, α = 0.1 / (2π), ΔE = 0.5)
end

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
	jldsave("benchmark/Graphene/ResultsNonSmearing/ValuesLT_N$(tabN)_final.jld2"; LTvalues, LTtimes, Energies, tabN, exDOS)
end



PTRvalues = zeros(length(Energies), length(ηlist), length(tabN))
PTRtimes = zeros(length(Energies), length(ηlist), length(tabN))

for (ie, E) in enumerate(Energies)
	@time for (iN, N) in enumerate(tabN)
		for (iη, η) in enumerate(ηlist)
			temp = @timed dos_solver_ptr(N, η)(E)
			PTRvalues[ie, iη, iN] = -imag(temp.value) / π / det(bz.B)
			PTRtimes[ie, iη, iN] = temp.time
		end
	end
	jldsave("benchmark/Graphene/ResultsNonSmearing/ValuesPTR_N$(tabN)_final2.jld2"; PTRvalues, PTRtimes, Energies, tabN, ηlist)
end

IAIvalues = zeros(length(Energies), length(ηlist), length(tolerances))
IAINs = zeros(length(Energies), length(ηlist), length(tolerances))
IAItimes = zeros(length(Energies), length(ηlist), length(tolerances))
for (itol, tol) in enumerate(tolerances)
	@time for (iη, η) in enumerate(ηlist)
		for (iE, E) in enumerate(Energies)
			temp = @timed dos_solver_iai(η, tol)(E)
			IAIvalues[iE, iη, itol] = -imag(temp.value[1]) / π / det(bz.B)
			IAINs[iE, iη, itol] = temp.value[2]
			IAItimes[iE, iη, itol] = temp.time
		end
	end
	jldsave("benchmark/Graphene/ResultsNonSmearing/ValuesIAI_final.jld2"; IAIvalues, IAINs, IAItimes, Energies, ηlist, tolerances)
end

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
jldsave("benchmark/Graphene/ResultsNonSmearing/ErrorsIAI_final.jld2"; IAIvalues, IAINs, IAItime, Energies, ηlist, tolerances, allerrorsIAI, tabNIAI, alltabtol, alltabη)

