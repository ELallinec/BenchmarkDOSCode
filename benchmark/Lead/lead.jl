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

using OffsetArrays
nb = 3

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

Energies = range(-1.0, 20.0, 301)

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
jldsave("benchmark/Lead/Results/BCD_N=100.0_E=-1.0-20.0.jld2"; BCDvalues, Energies, α, ΔE = 0.3)

LTvalues = zeros(length(Energies))
@time cache2 = AutoBZCore.init(probLT, AutoBZCore.LT(; npt = 100));
cache2.domain = Energies[50]
@time tempLT = AutoBZCore.solve!(cache2).value
@time for (ie, e) in enumerate(Energies[1:end])
	cache2.domain = e
	LTvalues[ie] = AutoBZCore.solve!(cache2).value
end
jldsave("benchmark/Lead/Results/LT_N=100.0_E=-1.0-20.0_full.jld2"; LTvalues, Energies)
