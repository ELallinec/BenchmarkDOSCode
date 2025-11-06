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
const J = 8;

using WannierIO

hrdat = read_w90_hrdat("benchmark/Silicon/silicon_hr.dat");
eigdat = read_w90_band_dat("benchmark/Silicon/silicon_band.dat");
eigkpt = read_w90_band_kpt("benchmark/Silicon/silicon_band.kpt");
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
ax.plot(bands[7, :], linewidth = 3, color = soft_olive)
ax.plot(bands[8, :], linewidth = 3, color = clear_gray)
ax.set_xlabel("k-points")
ax.set_xticks([1, 101, 216, 258, 380])
ax.set_xticklabels(["L", "G", "X", "K", "G"])
ax.vlines([1, 101, 216, 258, 380], -6, 17, colors = :black)

trueeig = reduce(hcat, eigdat.eigenvalues)
ax.plot(trueeig[1, :], linewidth = 3, color = gold, linestyle = :dashed)
ax.plot(trueeig[2, :], linewidth = 3, color = dark_navy, linestyle = :dashed)
ax.plot(trueeig[3, :], linewidth = 3, color = rust_orange, linestyle = :dashed)
ax.plot(trueeig[4, :], linewidth = 3, color = teal, linestyle = :dashed)
ax.plot(trueeig[5, :], linewidth = 3, color = slate_gray, linestyle = :dashed)
ax.plot(trueeig[6, :], linewidth = 3, color = light_blue, linestyle = :dashed)
ax.plot(trueeig[7, :], linewidth = 3, color = soft_olive, linestyle = :dashed)
ax.plot(trueeig[8, :], linewidth = 3, color = clear_gray, linestyle = :dashed)



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



Energies = range(-5.9, 16.5, 200)
bz = load_bz(FBZ(), I(d))
prob = DOSProblem(H, 0.99, bz)

tempEnergies = range(12, 13, 51)
#tempEnergies = Energies



probLT = DOSProblem(H1, 0.99, bz)

p0 = (; η = 1e-1, ω = 0.01)
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = AutoBZCore.FourierSeriesEvaluators.period(H1)
	greens_function(k, H1(k), p0)
end

integrand = FourierIntegralFunction(greens_function, H1, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, load_bz(CubicSymIBZ(), I(d)), p0; abstol = 1e-3)


function dos_solver_iai(η, atol)
	solver = init(AutoBZProblem(TrivialRep(), integrand, load_bz(CubicSymIBZ(), I(d)), p0; abstol = atol), EvalCounter(IAI()))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		@SVector [-imag(getproperty(temp, :value)) / det(bz.B) / π, getproperty(getproperty(temp, :stats), :numevals)]
	end
end

function dos_solver_ptr(N, η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, load_bz(CubicSymIBZ(), I(d)), p0), PTR(npt = N))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		-imag(getproperty(temp, :value) / det(bz.B)) / π
	end
end

IAIvalues = []
@time for E in tempEnergies
	temp = getindex(dos_solver_iai(1e-1, 1e-2)(E), 1)
	push!(IAIvalues, temp)
end
LTvalues = zeros(length(tempEnergies))
@time cache2 = AutoBZCore.init(prob, AutoBZCore.LT(; npt = 200));
cache2.domain = tempEnergies[25]
@time tempLT = AutoBZCore.solve!(cache2).value
@time for (ie, e) in enumerate(tempEnergies[1:end])
	cache2.domain = e
	LTvalues[ie] = AutoBZCore.solve!(cache2).value
end
jldsave("benchmark/Silicon/SiliconLT100E12-13.jld2"; Energies = tempEnergies, LTvalues)

BCDvalues = zeros(ComplexF64, length(tempEnergies))
@time cache1 = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = 30, α = 0.04 / 2π, ΔE = 0.3));
cache1.domain = tempEnergies[25]
@time tempBCD = AutoBZCore.solve!(cache1).value
@time for (ie, e) in enumerate(tempEnergies[1:end])
	cache1.domain = e
	BCDvalues[ie] = AutoBZCore.solve!(cache1).value
end
jldsave("benchmark/Silicon/SiliconBCD20E12-13.jld2"; Energies = tempEnergies, BCDvalues)