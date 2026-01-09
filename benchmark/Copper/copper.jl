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

hrdat = read_w90_hrdat("benchmark/Copper/Copper4x4x4/copper_hr.dat");
eigdat = read_w90_band_dat("benchmark/Copper/Copper4x4x4/copper_band.dat");
eigkpt = read_w90_band_kpt("benchmark/Copper/Copper4x4x4/copper_band.kpt");

using OffsetArrays
nb = 3

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

H = FourierSeries(H_R, period = 1)
H1 = HermitianFourierSeries(H)
DH = HessianSeries(H)
DH1 = HessianSeries(H1)
G = [0.000, 0.000, 0.000]
X = [0.500, 0.500, 0.000]
W = [0.500, 0.750, 0.250]
L = [0.000, 0.500, 0.000]
K = [0.000, 0.500, -0.500]

plot(reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))', lw = 3)
plot(reduce(hcat, eigdat.eigenvalues)', "--", lw = 3)
xticks([1, 101, 151, 222, 309, 450], ["G", "X", "W", "L", "G", "K"])
vlines([1, 101, 151, 222, 309, 450], minimum(reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))) * 0.99, maximum(reduce(hcat, eigdat.eigenvalues)) * 1.01, color = :k, lw = 3)
ylim(minimum(reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))) * 0.99, maximum(reduce(hcat, eigdat.eigenvalues)) * 1.01)
ylim(9.6, 10)
Energies = range(8, 9, 101)

bzBCD = load_bz(FBZ(), I(d))
bzLT = load_bz(FBZ(), I(d))
probBCD = DOSProblem(H, Energies[1], bzBCD)
probLT = DOSProblem(H, Energies[1], bzLT)

smallbz = [[x, y, z] for x in range(0, 1, 11), y in range(0, 1, 11), z in range(0, 1, 11)]
@time α = 2 / maximum(norm.(DH1.(smallbz)))
α = 0.003

BCDvalues = zeros(length(Energies))
@time cache1 = AutoBZCore.init(probBCD, AutoBZCore.BCD(; npt = 30, α = α, ΔE = 0.1));
cache1.domain = Energies[50]
@time tempBCD = AutoBZCore.solve!(cache1).value

@time for (ie, e) in enumerate(Energies[1:end])
	cache1.domain = e
	BCDvalues[ie] = AutoBZCore.solve!(cache1).value
end


LTvalues = zeros(length(Energies))
@time cache2 = AutoBZCore.init(probLT, AutoBZCore.LT(; npt = 50));
cache2.domain = Energies[50]
@time tempLT = AutoBZCore.solve!(cache2).value
@time for (ie, e) in enumerate(Energies[1:end])
	cache2.domain = e
	LTvalues[ie] = AutoBZCore.solve!(cache2).value
end
bands = reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))
gold = "#ffa600"
dark_navy = "#003f5c"
rust_orange = "#bc5090"
teal = "#2A9D8F"
slate_gray = "#7a5195"
light_blue = "#4393c3"
fig, ax = subplots(1, 2)
fig.subplots_adjust(wspace = 0.0)
ax[1].plot(bands[1, :], linewidth = 3, color = gold)
ax[1].plot(bands[2, :], linewidth = 3, color = dark_navy)
ax[1].plot(bands[3, :], linewidth = 3, color = rust_orange)
ax[1].plot(bands[4, :], linewidth = 3, color = teal)
ax[1].plot(bands[5, :], linewidth = 3, color = slate_gray)
ax[1].plot(bands[6, :], linewidth = 3, color = light_blue)

ax[1].set_xlabel("k-points")
ax[1].set_xticks([1, 101, 151, 222, 309, 450], ["G", "X", "W", "L", "G", "K"])
ax[1].vlines([1, 101, 151, 222, 309, 450], minimum(reduce(hcat, eigvals.(H1.(eigkpt.kpoints)))) * 0.99, maximum(reduce(hcat, eigdat.eigenvalues)) * 1.01, color = :k, lw = 3)
#ax[1].hlines(-8.787, 1, 1501, color = :b, linewidth = 2, label = "van Hove", linestyle = :dashed)
ax[1].set_ylim(8, 9.0)
ax[1].set_ylabel("Energies (eV)")


ax[2].plot(BCDvalues, Energies, color = :red, linewidth = 3, label = "BCD")

ax[2].plot(LTvalues, Energies, color = :black, linewidth = 3, label = "Reference", linestyle = :dotted)
ax[2].set_xticks([])
ax[2].set_xlabel("DOS")
ax[2].set_yticks([])
ax[2].set_xlim(0, 2.0)
ax[2].set_ylim(8, 9)
#ax[2].legend()
ax[2].fill_betweenx(Energies, BCDvalues, color = :black, alpha = 0.2)
ax[2].fill_betweenx(Energies, LTvalues, color = :black, alpha = 0.1)

tempBCDvalues = copy(BCDvalues)
@load "benchmark/Copper/CopperBCD100DE=0.1a=0.01E9.6-10.jld2"
