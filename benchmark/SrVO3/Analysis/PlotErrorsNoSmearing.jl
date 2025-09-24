using Pkg
Pkg.activate(".")
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 32)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using JLD2
path = "benchmark/SrVO3/ResultsNonSmearing/"
d = 3
exDOS, Energies = load(path * "exDOS_N=501.jld2", "exDOS", "Energies")
BCDvalues, BCDtimes = load(path * "ValuesBCD_N=3:200_final.jld2", "BCDvalues", "BCDtimes")
LTvalues, LTtimes = load(path * "ValuesLT_N=3:200_final.jld2", "LTvalues", "LTtimes")
IAIvalues, IAItimes, IAINs, tolerances, ηlist = load(path * "ValuesIAI_final2.jld2", "IAIvalues", "IAItimes", "IAINs", "tolerances", "ηlist")
PTRvalues, PTRtimes = load(path * "ValuesPTR_N=3:400_final.jld2", "PTRvalues", "PTRtimes")

PTRindices = [[argmin(abs.(PTRvalues[i, :, j] .- exDOS[i])) for j in eachindex(3:400)] for i in eachindex(Energies)]
PTRvalue = [PTRvalues[i, PTRindices[i][j], j] for i in eachindex(Energies), j in eachindex(3:400)]
PTRtime = [PTRtimes[i, PTRindices[i][j], j] for i in eachindex(Energies), j in eachindex(3:400)]
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

for i in eachindex(Energies)
	fig, ax = subplots()
	ax.scatter(PTRtime[i, 1:2:end], abs.(PTRvalue[i, 1:2:end] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "+", color = (0, 158 / 255, 115 / 255), label = "PTR", s = 150)
	ax.scatter(IAItime[i][1:2:end], abs.(allerrorsIAI[i][1:2:end]) ./ (abs(exDOS[i])), marker = "^", color = :black, label = "IAI", s = 100)
	ax.scatter(LTtimes[i, 1:2:end], abs.(LTvalues[i, 1:2:end] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "2", color = (230 / 255, 159 / 255, 0), label = "LT", s = 150)
	ax.scatter(BCDtimes[i, 1:2:end], abs.(BCDvalues[i, 1:2:end] .- exDOS[i]) ./ (abs(exDOS[i])), marker = "v", color = (204 / 255, 121 / 255, 167 / 255), label = "BCD", s = 100)
	
	ax.set_yscale(:log)
	ax.set_xscale(:log)
	#ax.set_ylim(1e-6,1)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Relative error")
	ax.set_title("E=$(Energies[i])")
	ax.legend()
	#savefig("benchmark/SrVO3/fig/SrVO3_Errors_vs_Time_E=$(Energies[i]).pdf")
end

close("all")

#=
f(k) = cos(k)
f_pt(N, M) = reduce(vcat, [cos(2π * k / N) * ones(M) for k in 0:N-1])
tabk = range(0, 2π, 1000)
plot(tabk, f.(tabk), "-k", linewidth = 2)
plot(tabk, f_pt(10, 100), "-r", linewidth = 2)
xlabel("k")

tab = range(-π, π, 11)[1:end-1]
fig, ax = subplots()
for x in tab
	ax.scatter(tab, x * ones(length(tab)), s = 80, color = :black)
end

ax.vlines(-π, -π, π, color = :red, linewidth = 2)
ax.vlines(π, -π, π, color = :red, linewidth = 2)
ax.hlines(-π, -π, π, color = :red, linewidth = 2)
ax.hlines(π, -π, π, color = :red, linewidth = 2)
ax.set_xlabel(L"k_x")
ax.set_ylabel(L"k_y")
ax.set_title(L"\text{Discretized Brillouin zone } \mathcal{B}_N")


e(k) = k^2
def(k, E) = -ForwardDiff.derivative(e, k) * exp(-(e(k) - E) / 0.5^2)
E = 1.5
η = 1.0
pole = sqrt(E + im * η)
contour(θ) = imag(pole) * sin(θ)
cont = vcat(zeros(136), -contour.(range(π, 2π, 318)), zeros(92), -contour.(range(0, π, 318)), zeros(136))
fig, ax = subplots()
ax.plot(tabk .- π, zeros(1000), "-k", linewidth = 1)
ax.plot(tabk .- π, cont, "--r", linewidth = 3)

ax.scatter(real(pole), imag(pole), color = :black, marker = :o, s = 80)
ax.scatter(-real(pole), -imag(pole), color = :black, marker = :o, s = 80)

ax.text(real(pole) + 0.1, imag(pole) + 0.1, L"\sqrt{E+i\eta}")
ax.text(-real(pole) - 1.3, -imag(pole) - 0.2, L"-\sqrt{E+i\eta}")
ax.arrow(-real(pole) + 0.2, -imag(pole) - 0.1, 0, 0.2, head_width = 0.05, color = :black)
ax.arrow(real(pole) - 0.2, imag(pole) + 0.1, 0, -0.2, head_width = 0.05, color = :black)
#legend()
ax.set_ylim(-1.2, 1.2)
=#