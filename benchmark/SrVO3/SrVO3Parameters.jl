
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
const J = 3;

using WannierIO

hrdat = read_w90_hrdat("data/svo_hr.dat");

Rmin, Rmax = extrema(hrdat.Rvectors)
Rsize = Tuple(Rmax .- Rmin .+ 1)
n, m = size(first(hrdat.H))

using OffsetArrays

H_R = OffsetArray(
	Array{SMatrix{n, m, eltype(eltype(hrdat.H)), n * m}}(undef, Rsize...),
	map(:, Rmin, Rmax)...,
)

allR = OffsetArray(
	Array{SVector{3}}(undef, Rsize...),
	map(:, Rmin, Rmax)...,
)
for (i, h, n) in zip(hrdat.Rvectors, hrdat.H, hrdat.Rdegens)
	H_R[CartesianIndex(Tuple(i))] = h / n
	allR[CartesianIndex(Tuple(i))] = [Tuple(i)[1], Tuple(i)[2], Tuple(i)[3]]
end



H = AutoBZCore.FourierSeries(H_R, period = 1);

dH = AutoBZCore.HessianSeries(H)

ηlist = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]), rev = true)               # 10 meV (scattering amplitude)
tolerances = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:4]), rev = true)
Energies = range(11, 14, 101)

p0 = (; η = 1e-2, ω = 12.44) # initial parameters
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end
integrand = FourierIntegralFunction(greens_function, H, prototype)

function dos_solver_ptr(N, η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz1, p0), PTR(npt = N))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end