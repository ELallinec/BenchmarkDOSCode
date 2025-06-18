
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

@load "benchmark/SrVO3/Results/RefDos.jld2"

#temp=map(i -> maximum([norm(BCDvalues[i, 1:77] - RefDos[1:77], Inf), norm(BCDvalues[i, 79:end] - RefDos[79:end], Inf)]), 1:63)
