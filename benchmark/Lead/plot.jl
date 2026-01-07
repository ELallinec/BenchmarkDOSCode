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
nb = 2

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

N=100

bands = reduce(hcat, eigdat.eigenvalues)
wbands=reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))
gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"

fig, ax = subplots()
ax.plot(bands[1, :], linewidth = 3, color = gold)
ax.plot(bands[2, :], linewidth = 3, color = dark_navy)
ax.plot(bands[3, :], linewidth = 3, color = rust_orange)
ax.plot(bands[4, :], linewidth = 3, color = teal)

ax.plot(wbands[1, :], linewidth = 3, linestyle = :dashed)
ax.plot(wbands[2, :], linewidth = 3, linestyle = :dashed)
ax.plot(wbands[3, :], linewidth = 3, linestyle = :dashed)
ax.plot(wbands[4, :], linewidth = 3, linestyle = :dashed)
ax.set_xlabel("k-points")
#ax.set_xticks([1, 3*N, 7*N, 9*N, 12*N], ["L", "G", "X", "K", "G"])
#ax.vlines([1, 3*N, 7*N, 9*N, 12*N], minimum(bands) * 0.99, maximum(bands) * 1.01, color = :k, lw = 3)


fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], linewidth = 3, color = gold)
ax[1].plot(bands[2, :], linewidth = 3, color = dark_navy)
ax[1].plot(bands[3, :], linewidth = 3, color = rust_orange)
ax[1].plot(bands[4, :], linewidth = 3, color = teal)


ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 101, 151, 222, 309, 450], ["G", "X", "W", "L", "G", "K"])
#ax[1].vlines([1, 101, 151, 222, 309, 450], minimum(reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))) * 0.99, maximum(reduce(hcat, eigdat.eigenvalues)) * 1.01, color = :k, lw = 3)
#ax[1].hlines(-8.787, 1, 1501, color = :b, linewidth = 2, label = "van Hove", linestyle = :dashed)
#ax[1].set_ylim(8, 9.0)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(BCDvalues, Energies, color = :red, linewidth = 3, label = "BCD")

ax[2].plot(LTvalues, Energies, color = :black, linewidth = 3, label = "Reference", linestyle = :dotted)
ax[2].set_xticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(0, 2.0)
#ax[2].set_ylim(8, 9)
#ax[2].legend()
ax[2].fill_betweenx(Energies, BCDvalues, color = :black, alpha = 0.2)
ax[2].fill_betweenx(Energies, LTvalues, color = :black, alpha = 0.1)
