using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 28)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2

@load "benchmark/Mono1D/Results/ValuesBCD_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesPTR_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesLT_N3:1500_final.jld2"
@load "benchmark/Mono1D/Results/ValuesIAI_final.jld2"

#BCDtime = [minimum(BCDtimes[:, iN]) for iN in 1:1498]
PTRtime = [[minimum(PTRtimes[i, :, iN]) for iN in 1:1498] for i in eachindex(Energies)]
PTRerror = [minimum(abs.(PTRvalues[i, :, iN] .- exDOS[i])) for i in 1:5, iN in 1:1498]


#! IAI absolute error
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
#jldsave(repo * "ErrorsIAI.jld2"; IAIvalues, IAINs, Energies, ηlist, tolerances, allerrorsIAI, tabNIAI, alltabtol, alltabη, exDOS)


for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(PTRtime[i, :], PTRerror[i, :] ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR", s = 150)
	ax.scatter(tempIAItimes[i,], allerrorsIAI[i] ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI", s = 150)
	ax.scatter(LTtimes[i, :], abs.(LTvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT", s = 150)
	ax.scatter(BCDtimes[i, :], abs.(BCDvalues[i, :] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", s = 150)
	yscale(:log)
	xlim(0, 1e-4)
	xlabel("Time (s)")
	ylabel("Relative error")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	#savefig("benchmark/Mono1D/fig/Mono1D_RelErrorsVsTime_E=$(Energies[i]).pdf")
end
close("all")