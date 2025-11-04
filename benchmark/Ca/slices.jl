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
using Revise
#Common parameters
const d = 3;
const J = 4;

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
# Example: kz = 0 plane
using LinearAlgebra, GLMakie, Base.Threads


E = (x, y, z) -> real(eigvals(H1([x, y, z])))

p = [0, 0, 0.0]
no = [1.0, 1.0, 0.0]
function compute_crossing(p, no, N = 101)
	M = nullspace(reshape(no, 1, 3))
	a, b = M[:, 1], M[:, 2]

	ul = vl = range(-0.5, 0.5, length = N)[1:end-1]
	U = [ul[j] for i in 1:(N-1), j in 1:(N-1)]
	V = [vl[i] for i in 1:(N-1), j in 1:(N-1)]
	X = @. a[1] * U + b[1] * V + p[1]
	Y = @. a[2] * U + b[2] * V + p[2]
	Z = @. a[3] * U + b[3] * V + p[3]

	# wrap coordinates into [0,1] using modulus for periodicity
	X = X
	Y = Y
	Z = Z

	surfcolor = Array{Float32}(undef, N - 1, N - 1)
	@threads for i in 1:N-1
		for j in 1:N-1
			vals = E(X[i, j], Y[i, j], Z[i, j])
			surfcolor[i, j] = Float32(minimum(abs.(diff(sort(vals)))))
		end
	end



	return surfcolor, X, Y, Z
end

function compute_crossings(surfcolor, tol = 5e-5)
	xs = Float32[]
	ys = Float32[]
	zs = Float32[]
	for i in axes(surfcolor, 1), j in axes(surfcolor, 2)
		if surfcolor[i, j] < tol
			push!(xs, X[i, j])
			push!(ys, Y[i, j])
			push!(zs, Z[i, j])
		end
	end
	return xs, ys, zs
end
G = zeros(3)
M = 0.5 * [1, 0, 0]
K = (1 / 3) * [1, 1, 0]
A = [0, 0, 0.5]
L = [0.5, 0, 0.5]
Ha = [1 / 3, 1 / 3, 1 / 2]
p = K
no = [1, 1, 0]
@time surfcolor, X, Y, Z = compute_crossing(p, no, 101)
@time xs, ys, zs = compute_crossings(surfcolor, 1e-3)
fig = Figure(resolution = (1800, 1300))
ax = Axis3(fig[1, 1])
surface!(ax, X, Y, Z, color = fill(RGBA(1, 1, 1, 0.1), size(Z)))
scatter!(ax, xs, ys, zs, markersize = 4, color = :black)
fig



temp = zeros(SVector{4}, (100, 100))
@time for ix in eachindex(X)
	temp[ix] = real(eigvals(H1([X[ix], Y[ix], Z[ix]])))
end

u = [1, -1, 0]
v = [0, 0, 1]
temp = zeros(SVector{4}, 100, 100)
Xs = range(0, 1, 101)[1:end-1]
for (ix, tx) in enumerate(range(0, 1, 101)[1:end-1])
	for (iy, ty) in enumerate(range(0, 1, 101)[1:end-1])

		temp[ix, iy] = real(eigvals(H1(tx * u + ty * v + K)))
	end
end
fig = Figure(resolution = (1800, 1300))
ax = Axis3(fig[1, 1])
surface!(ax, Xs, Xs, getindex.(temp, 1), color = fill(RGBA(1, 0, 0, 0.5), size(getindex.(temp, 1))))
surface!(ax, Xs, Xs, getindex.(temp, 2), color = fill(RGBA(0, 1, 0, 0.5), size(getindex.(temp, 1))))
surface!(ax, Xs, Xs, getindex.(temp, 3), color = fill(RGBA(0, 0, 1, 0.5), size(getindex.(temp, 1))))
surface!(ax, Xs, Xs, getindex.(temp, 4), color = fill(RGBA(1, 0, 1, 0.5), size(getindex.(temp, 1))))
fig