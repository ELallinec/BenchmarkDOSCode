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

#Common parameters
const d = 3;
const J = 5;

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
#=
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
ylim(6.5, 8.25)

=#
Energies = range(9.6, 10, 101)

bzBCD = load_bz(FBZ(), I(d))
bzLT = load_bz(FBZ(), I(d))
probBCD = DOSProblem(H, Energies[1], bzBCD)
probLT = DOSProblem(H, Energies[1], bzLT)
#=
α = 0.03 / 2π
ΔE = 0.05

BCDvalues = zeros(length(Energies))
@time cache1 = AutoBZCore.init(probBCD, AutoBZCore.BCD(; npt = 3, α = α, ΔE = ΔE));
cache1.domain = Energies[50]
@time tempBCD = AutoBZCore.solve!(cache1).value
@time cache1 = AutoBZCore.init(probBCD, AutoBZCore.BCD(; npt = 100, α = α, ΔE = ΔE));
@time for (ie, e) in enumerate(Energies[1:end])
	cache1.domain = e
	BCDvalues[ie] = AutoBZCore.solve!(cache1).value
end

jldsave("benchmark/Copper/CopperBCD100DE=$(ΔE)a=$(α*2π)E9.6-10.jld2"; BCDvalues, Energies, α, ΔE)
=#
LTvalues = zeros(length(Energies))
@time cache2 = AutoBZCore.init(probLT, AutoBZCore.LT(; npt = 3));
cache2.domain = Energies[50]
@time tempLT = AutoBZCore.solve!(cache2).value
@time cache2 = AutoBZCore.init(probLT, AutoBZCore.LT(; npt = 500));
@time for (ie, e) in enumerate(Energies[1:end])
	cache2.domain = e
	LTvalues[ie] = AutoBZCore.solve!(cache2).value
end

jldsave("benchmark/Copper/CopperLT700.jld2"; LTvalues, Energies)
