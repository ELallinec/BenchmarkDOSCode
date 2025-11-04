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
using Colors
#using DensityOfStates
using ForwardDiff
using FourierSeriesEvaluators
using Revise
#Common parameters
const d = 3;
const J = 4;
using PyPlot
PyPlot.rc("font", family = "serif")
PyPlot.rc("mathtext", fontset = "dejavuserif")
PyPlot.rc("font", size = 25)
PyPlot.rc("figure", figsize = (9 * 1.5, 6 * 9 * 1.5 / 8))
using WannierIO

hrdat = read_w90_hrdat("benchmark/Ca/OldResults/Ca_hr.dat");
using OffsetArrays
nb = 8
#bandss = WannierIO.read_qe_band("benchmark/Ca/Ca_bands.dat")
Rvecs = hrdat.Rvectors

tempR = []
for (iR, R) in enumerate(Rvecs)
	if abs(R[1]) <= nb && abs(R[2]) <= nb && abs(R[3]) <= nb
		push!(tempR, iR)
	end
end

H_R = OffsetArray(zeros(SMatrix{J, J, ComplexF64}, 2 * nb + 1, 2 * nb + 1, 2 * nb + 1), -nb:nb, -nb:nb, -nb:nb)
C = hrdat.H
for i in tempR
	H_R[Rvecs[i][1], Rvecs[i][2], Rvecs[i][3]] = C[i] / hrdat.Rdegens[i]
end
H1 = HermitianFourierSeries(FourierSeries(H_R, period = 1))
H = FourierSeries(H_R, period = 1)

G = zeros(3)
M = 0.5 * [1, 0, 0]
K = (1 / 3) * [1, 1, 0]
A = [0, 0, 0.5]
L = [0.5, 0, 0.5]
Ha = [1 / 3, 1 / 3, 1 / 2]
tabpath = [G, M, K, G, A, L, Ha]
stringpath = ["G", "M", "K", "G", "A", "L", "H"]
function segment(p1, p2, N)
	return [(1 - n / N) * p1 + n / N * p2 for n in 0:N]
end

function path(tabpath, N)
	return reduce(vcat, [(i == length(tabpath) - 1 ? segment(tabpath[i], tabpath[i+1], N) : segment(tabpath[i], tabpath[i+1], N)[1:end-1]) for i in eachindex(tabpath[1:end-1])])
end
npath = 10000
thepath = path(tabpath, npath)
kpath_with_weights = [vcat((k), 1.0) for k in thepath]   # each k-point becomes [kx, ky, kz, 1.0]

@time bands = reduce(hcat, map(k -> real(eigen(H1((k))).values), thepath))
@load "benchmark/Ca/BCD50_all2.jld2"
@load "benchmark/Ca/LT300_all.jld2"
gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
soft_olive = "#a6a631"
clear_gray = "#c0c0c0"


fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], linewidth = 3, color = gold)
ax[1].plot(bands[2, :], linewidth = 3, color = dark_navy)
ax[1].plot(bands[3, :], linewidth = 3, color = rust_orange)
ax[1].plot(bands[4, :], linewidth = 3, color = teal)

ax[1].set_xlabel("k-points")
ax[1].set_xticks(npath * (0:6) .+ 1)
ax[1].set_xticklabels(["G", "M", "K", "G", "A", "L", "H"])
ax[1].vlines(npath * (0:6) .+ 1, 0, 8, colors = :black)
ax[1].hlines([6.56, 6.88], 1, npath * 6 + 1, colors = [:k, :k], linestyles = :dashed, linewidths = 3, label = "Crossing")
ax[1].hlines([5.175, 5.55], 1, npath * 6 + 1, colors = [:k, :k], linestyles = :dashed, linewidths = 3)
ax[1].set_ylim(0.9, 7.3)
ax[1].set_xlim(0, 6 * npath + 1)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(BCDvalues, Energies, color = :red, linewidth = 3, label = "BCD")
ax[2].plot(LTvalues, Energies, color = :black, linewidth = 3, label = "Reference", linestyle = :dotted)
ax[2].hlines([6.56, 6.88], 0, 2.5, colors = [:k, :k], linestyles = :dashed, linewidths = 3)
ax[2].hlines([5.175, 5.55], 0, 2.5, colors = [:k, :k], linestyles = :dashed, linewidths = 3)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_xlim(0, 2.3)
ax[2].set_ylim(0.9, 7.3)
#ax[2].legend()
ax[2].fill_betweenx(Energies, BCDvalues, color = :black, alpha = 0.2)
ax[2].fill_betweenx(Energies, LTvalues, color = :black, alpha = 0.1)
fig.legend(loc = "lower right", bbox_to_anchor = (0.87, 0.18))
#=
for i in 1:4

	ax[3].scatter(gaps[i, :], tab[i, :], s = 5, color = tabcolors[i, :])
end
ax[3].set_xlabel("Minimal gap")
ax[3].set_yticks([])
ax[3].set_ylim(0.9, 7.3)
=#


bz = load_bz(FBZ(), I(d))
probLT = prob = DOSProblem(H, 0.99, bz)
#tempEnergies = range(5, 7.3, 50)
tempEnergies = Energies
temps = []
for t in tabcolors[1, :]
	if t == :black
		push!(temps, 0)
	else
		push!(temps, 10)
	end
end



LTvalues = zeros(length(tempEnergies))
LTtimes = zeros(length(tempEnergies))
@time cache2 = AutoBZCore.init(prob, AutoBZCore.LT(; npt = 300));
@time for (ie, e) in enumerate(tempEnergies)
	cache2.domain = e
	temp = @timed AutoBZCore.solve!(cache2).value
	LTvalues[ie] = temp.value
	LTtimes[ie] = temp.time
end
jldsave("benchmark/Ca/LT300_all.jld2"; LTvalues, Energies = tempEnergies)
BCDvalues = zeros(length(tempEnergies))
BCDtimes = zeros(length(tempEnergies))
@time cache1 = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 100, ΔE = 0.2));
@time for (ie, e) in enumerate(tempEnergies[73:73])
	cache1.domain = e
	BCDvalues[ie] = AutoBZCore.solve!(cache1).value
end
jldsave("benchmark/Ca/BCD50_all2.jld2"; BCDvalues, Energies = tempEnergies)
function gap_intelligent(H, BZ, allgap)
	Abuf = @MMatrix zeros(ComplexF64, J, J)
	ev = @MVector zeros(Float64, J)
	temp = @MVector zeros(Float64, J - 1)
	# Compute all k-point gaps

	@inbounds for (ik, k) in enumerate(BZ)
		Abuf .= H(k)
		ev .= eigvals(SMatrix(Abuf))
		temp .= diff(ev)
		for j in 2:J-1
			allgap[ik][j] = min(temp[j-1], temp[j])
		end
		allgap[ik][1] = temp[1]
		allgap[ik][end] = temp[end]
	end
	#return allgap
end

bz_iter = Iterators.product(range(0.0, stop = 1.0, length = 26)[1:end-1], range(0.0, stop = 1.0, length = 26)[1:end-1], range(0.0, stop = 1.0, length = 26)[1:end-1])
smallbz = [SVector{3, Float64}(kx, ky, kz) for (kx, ky, kz) in bz_iter]
otherbz = [smallbz[ik] for ik in 1:length(smallbz)]

allgap = [@MVector zeros(Float64, J) for _ in 1:length(smallbz)]
@time gap_intelligent(H1, smallbz, allgap)
gaps = reduce(hcat, allgap)
@time tabeigen = eigen.(H1.(otherbz))
tabλ = getproperty.(tabeigen, :values)
tab = reduce(hcat, tabλ)
tabψ = getproperty.(tabeigen, :vectors)
tabh = H1.(otherbz)
@time tabdh = getindex.(DH.(otherbz), 2)
tabtr = [[tr(tabdh[ih][j] * exp(-(tabh[ih] - 4.5 * I)^2 / 0.2^2)) for j in 1:3] for ih in eachindex(tabdh)]
cbz = otherbz - im * 0.1 / 2π * tabtr
tabcλ = -imag([eigen(H(k)).values for k in cbz])
tabc = reduce(hcat, tabcλ)
tabcolors = Array{Symbol}(undef, size(tabc))
for index in CartesianIndices(tabc)
	if tabc[index] < -1e-10

		tabcolors[index] = colors5[index[1]+1]
	else
		tabcolors[index] = colors5[1]
	end
end


colors5 = [
	:black,    # 0
	:blue,      # 1
	:orange,     # 2
	:green,    # 3
	:red,   # 4
	# 15
]
uniquecomb = sort(unique(unique.(eachcol(tabcolors))))

temp = reduce(vcat, [colors5[findall(k -> (unique(tabcolors[:, ik]) == k || [unique(tabcolors[:, ik])] == k), uniquecomb)] for ik in eachindex(otherbz)])

