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
const J = 6;

using WannierIO

hrdat = read_w90_hrdat("benchmark/Silicon/silicon_hr.dat");
eigdat = read_w90_band_dat("benchmark/Silicon/silicon_band.dat");
eigkpt = read_w90_band_kpt("benchmark/Silicon/silicon_band.kpt");
using OffsetArrays
nb = 4

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

H = FourierSeries(H_R, period = 1)
H1 = HermitianFourierSeries(H)
DH = HessianSeries(H)
DH1 = HessianSeries(H1)
#! Verification de la tronche des bandes
L = [0.5, 0.5, 0.5]
G = [0.0, 0.0, 0.0]
X1 = [0.5, 0.0, 0.5]
X2 = [0.5, -0.5, 0.0]
K = [0.375, -0.375, 0.0]
tabpath = [L, G, X1, X2, K, G]

thepath = eigkpt.kpoints
bands = reduce(hcat, map(k -> eigvals(H1(k)), thepath))


gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
soft_olive = "#a6a631"
clear_gray = "#c0c0c0"

fig, ax = subplots()
fig.subplots_adjust(wspace = 0.0)
ax.plot(bands[1, :], linewidth = 3, color = gold)
ax.plot(bands[2, :], linewidth = 3, color = dark_navy)
ax.plot(bands[3, :], linewidth = 3, color = rust_orange)
ax.plot(bands[4, :], linewidth = 3, color = teal)
ax.plot(bands[5, :], linewidth = 3, color = slate_gray)
ax.plot(bands[6, :], linewidth = 3, color = light_blue)

ax.set_xlabel("k-points")
ax.set_xticks([1, 101, 216, 258, 380])
ax.set_xticklabels(["L", "G", "X", "K", "G"])
ax.vlines([1, 101, 216, 258, 380], minimum(trueeig)*0.99, maximum(trueeig)*1.01, colors = :black, lw = 3)

trueeig = reduce(hcat, eigdat.eigenvalues)
ax.plot(trueeig[1, :], linewidth = 3, color = gold, linestyle = :dashed)
ax.plot(trueeig[2, :], linewidth = 3, color = dark_navy, linestyle = :dashed)
ax.plot(trueeig[3, :], linewidth = 3, color = rust_orange, linestyle = :dashed)
ax.plot(trueeig[4, :], linewidth = 3, color = teal, linestyle = :dashed)
ax.plot(trueeig[5, :], linewidth = 3, color = slate_gray, linestyle = :dashed)
ax.plot(trueeig[6, :], linewidth = 3, color = light_blue, linestyle = :dashed)



fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], linewidth = 3, color = gold)
ax[1].plot(bands[2, :], linewidth = 3, color = dark_navy)
ax[1].plot(bands[3, :], linewidth = 3, color = rust_orange)
ax[1].plot(bands[4, :], linewidth = 3, color = teal)
ax[1].plot(bands[5, :], linewidth = 3, color = slate_gray)
ax[1].plot(bands[6, :], linewidth = 3, color = light_blue)
ax[1].plot(bands[7, :], linewidth = 3, color = soft_olive)
ax[1].plot(bands[8, :], linewidth = 3, color = clear_gray)
ax[1].set_xlabel("k-points")
ax[1].set_xticks(vcat(1, cumsum(npt) .+ 1))
ax[1].set_xticklabels(["L", "G", "X", "K", "G"])
ax[1].vlines(vcat(1, cumsum(npt) .+ 1), -6, 17, colors = :black)
#ax[1].hlines(-8.787, 1, 1501, color = :b, linewidth = 2, label = "van Hove", linestyle = :dashed)

ax[1].set_ylabel("Energies (eV)")


ax[2].plot(BCDvalues, tempEnergies, color = :red, linewidth = 3, label = "BCD")

ax[2].plot(LTvalues, tempEnergies, color = :black, linewidth = 3, label = "Reference", linestyle = :dotted)
ax[2].hlines([-8.74, -8.622, -8.18, -7.728], 0, 0.7, colors = [:k, :k, :k, :k], linestyles = :dashed, linewidths = 3)
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
#ax[2].set_xlim(0, 0.65)
ax[2].set_ylim(-9.0, -7)
#ax[2].legend()
ax[2].fill_betweenx(tempEnergies, BCDvalues, color = :black, alpha = 0.2)
ax[2].fill_betweenx(tempEnergies, LTvalues, color = :black, alpha = 0.1)
fig.legend(loc = "center right", bbox_to_anchor = (0.9, 0.77))




Energies = range(6.5, 18.5, 301)

bzBCD = load_bz(FBZ(), I(d))
bzLT = load_bz(FBZ(), I(d))
probBCD = DOSProblem(H, Energies[1], bzBCD)
probLT = DOSProblem(H, Energies[1], bzLT)

smallbz = [[x, y, z] for x in range(0, 1, 11), y in range(0, 1, 11), z in range(0, 1, 11)]
α = 1 / maximum(norm.(DH1.(smallbz)))



BCDvalues = zeros(length(Energies))

@time cache1 = AutoBZCore.init(probBCD, AutoBZCore.BCD(; npt = 30, α = α, ΔE = 0.3));

@time tempBCD = AutoBZCore.solve!(cache1).value

@time for (ie, e) in enumerate(Energies[1:end])
	cache1.domain = e
	BCDvalues[ie] = AutoBZCore.solve!(cache1).value
end


LTvalues = zeros(length(Energies))
@time cache2 = AutoBZCore.init(probLT, AutoBZCore.LT(; npt = 100));
cache2.domain = Energies[50]
@time tempLT = AutoBZCore.solve!(cache2).value
@time for (ie, e) in enumerate(Energies[1:end])
	cache2.domain = e
	LTvalues[ie] = AutoBZCore.solve!(cache2).value
end