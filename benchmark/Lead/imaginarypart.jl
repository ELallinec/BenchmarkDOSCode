using Pkg: Pkg
Pkg.activate(".")
push!(LOAD_PATH, "src")
#Pkg.instantiate()
using LinearAlgebra
using AutoBZCore
using FFTW
using Folds
using BenchmarkTools
using StaticArrays
using JLD2
#using DensityOfStates
using ForwardDiff
using FourierSeriesEvaluators
using PyPlot
using Colors
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
#Common parameters
const d = 3;
const J = 4;
using WannierIO

hrdat = read_w90_hrdat("benchmark/Lead/lead_hr.dat");
eigdat = read_w90_band_dat("benchmark/Lead/lead_band.dat");
eigkpt = read_w90_band_kpt("benchmark/Lead/lead_band.kpt");
using OffsetArrays
nb = 6

Rvecs = hrdat.Rvectors
tempR = []
for (iR, R) in enumerate(Rvecs)
	if abs(R[1]) <= nb && abs(R[2]) <= nb && abs(R[3]) <= nb
		push!(tempR, iR)
	end
end

H_R = OffsetArray(zeros(SMatrix{J, J, ComplexF64}, 2 * nb + 1, 2 * nb + 1, 2 * nb + 1), (-nb):nb, (-nb):nb, (-nb):nb)
C = hrdat.H
for i in tempR
	H_R[Rvecs[i][1], Rvecs[i][2], Rvecs[i][3]] = C[i] / hrdat.Rdegens[i]
end

h = FourierSeries(H_R, period = 1)
h1 = HermitianFourierSeries(h)
dh = HessianSeries(h)
dh1 = HessianSeries(h1)
bands = reduce(hcat, eigdat.eigenvalues)
wbands=reduce(hcat, eigvals.(h1.(eigkpt.kpoints)))
gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
black = "#000000"

bzBCD = load_bz(FBZ(), I(d))
wd2h = FourierSeriesEvaluators.workspace_allocate(dh, FourierSeriesEvaluators.period(dh))
kalg = MonkhorstPack(npt = 30, syms = bzBCD.syms)
dom = AutoBZCore.canonical_ptr_basis(bzBCD.B)
rule = AutoBZCore.init_fourier_rule(wd2h, dom, kalg)

if rule isa AutoBZCore.FourierPTR
	if typeof(rule.s[1][1]) <: Number
		tabeigen = eigen.(getindex.(rule.s, 1))
	else
		tabeigen = eigen.(Hermitian.(getindex.(rule.s, 1)))
	end
elseif rule isa FourierMonkhorstPack
	if typeof(rule.wxs[1][2].s) <: Number
		tabeigen = eigen.(getindex.(getproperty.(getindex.(rule.wxs, 2), :s), 1))
	else
		tabeigen = eigen.(Hermitian.(getindex.(getproperty.(getindex.(rule.wxs, 2), :s), 1)))
	end
end
α=0.001
η=0.0
E=9.57
ΔE = 0.05
tempdef=zeros(SVector{3, Float64}, size(rule.p))
tempddef = zeros(SMatrix{3, 3, Float64, 9}, size(rule.p))
@time for index in CartesianIndices(rule.p)
	H, d1H, d2H = rule.s[index]
	k = rule.p[index][2]
	tempdef[index], tempddef[index] = real.(AutoBZCore.return_def_and_ddef(h, k, H, d1H, d2H, tabeigen[index].values, tabeigen[index].vectors, Val(d), Val(J), E, α, η, ΔE))

end


function gap_intelligent(H, BZ, allgap)
	Abuf = @MMatrix zeros(ComplexF64, J, J)
	ev = @MVector zeros(Float64, J)
	temp = @MVector zeros(Float64, J - 1)
	# Compute all k-point gaps

	@inbounds for (ik, k) in enumerate(BZ)
		Abuf .= H(k)
		ev .= eigvals(SMatrix(Abuf))
		temp .= diff(ev)
		for j in 2:(J-1)
			allgap[ik][j] = min(temp[j-1], temp[j])
		end
		allgap[ik][1] = temp[1]
		allgap[ik][end] = temp[end]
	end
	#return allgap
end


smallbz = getindex.(rule.p, 2)
otherbz = [smallbz[ik] for ik in 1:length(smallbz)]

allgap = [@MVector zeros(Float64, J) for _ in 1:length(smallbz)]
@time gap_intelligent(h1, smallbz, allgap)
gaps = reduce(hcat, allgap)
#@time tabeigen = eigen.(h1.(otherbz))
tabλ = getproperty.(tabeigen, :values)
@time tab = reduce(hcat, [tabλ[ik] for ik in 1:length(tabλ)])
#tabψ = getproperty.(tabeigen, :vectors)
#tabh = h1.(otherbz)
#@time tabdh = getindex.(dh.(otherbz), 2)
#tabtr = [[tr(tabdh[ih][j] * exp(-(tabh[ih] - 9.5 * I)^2 / 0.2^2)) for j in 1:3] for ih in eachindex(tabdh)]
#cbz = otherbz - im * 0.1 / 2π * tabtr
cbz1 = smallbz .- α*im*tempdef
cbz = [cbz1[ik] for ik in 1:length(cbz1)]
@time tabcλ = -imag([eigen(h(k)).values for k in cbz])
tabc=reduce(hcat, tabcλ)
tabcolors = Array{String}(undef, size(tabc))
colors5 = [gold, dark_navy, rust_orange, teal]
for index in CartesianIndices(tabc)
	if tabc[index] < -1e-10
		tabcolors[index] = colors5[index[1]]

	else
		tabcolors[index] = black
	end
end


uniquecomb = sort(unique(unique.(eachcol(tabcolors))))


εf=11.2195


fig, ax = subplots(1, 2, width_ratios = (2, 1))
fig.subplots_adjust(wspace = 0.0)
ax[1].set_ylabel("Energy")
ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 101, 151, 222, 309, 450], ["G", "X", "W", "L", "G", "K"])
ax[1].vlines([1, 101, 151, 222, 309, 450], minimum(reduce(hcat, wbands) .- εf) * 1.05, maximum(reduce(hcat, wbands) .- εf) * 1.05, color = :k, lw = 3)
ax[2].set_xlabel("Gap")
ax[2].set_yticks([])
ax[1].set_ylim(minimum(reduce(hcat, wbands) .- εf) * 1.05, maximum(reduce(hcat, wbands) .- εf) * 1.05)
ax[2].set_ylim(minimum(reduce(hcat, wbands) .- εf) * 1.05, maximum(reduce(hcat, wbands) .- εf) * 1.05)

ax[1].plot(wbands[1, :] .- εf, color = gold, lw = 5)
ax[1].plot(wbands[2, :] .- εf, color = dark_navy, lw = 5)
ax[1].plot(wbands[3, :] .- εf, color = rust_orange, lw = 5)
ax[1].plot(wbands[4, :] .- εf, color = teal, lw = 5)

mask_black = tabcolors .== black
mask_other = .!mask_black

# Points noirs en arrière-plan
ax[2].scatter(gaps[mask_black], tab[mask_black] .- εf, c = "black", alpha = 0.9, label = "Non-negative", s = 30)

# Autres couleurs par-dessus
ax[2].scatter(gaps[mask_other], tab[mask_other] .- εf, c = tabcolors[mask_other], s = 30)
fig.legend(loc = "upper right", markerscale = 2.5)

using PyCall

@pyimport matplotlib.patches as mpatches
@pyimport matplotlib.transforms as mtransforms
function multicolor_legend_handle(colors; size = 10)
	fig = gcf()
	ax = gca()
	handles = PyObject[]

	for (i, c) in enumerate(colors)
		h = mlines.Line2D([], [], marker = "o", linestyle = "",
			markerfacecolor = c,
			markersize = size)
		push!(handles, h)
	end
	return tuple(handles...)
end

h_neg = multicolor_legend_handle((gold, dark_navy, rust_orange, teal))

fig.legend(
	Any[h0, h_neg],
	["Non-negative"],
	loc = "upper right",
)

fig.legend(
	Any[h1, h2, h3, h4],
	["Negative", "Negative", "Negative", "Negative"],
	loc = "upper right",
)
