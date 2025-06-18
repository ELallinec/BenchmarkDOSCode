using Pkg: Pkg

push!(LOAD_PATH, "src")
using LinearAlgebra
using AutoBZCore
using FFTW
using Folds
using BenchmarkTools
using StaticArrays
using JLD2
using DensityOfStates
using ForwardDiff
using OffsetArrays
using Elliptic
#Common parameters
d = 2
t = 1
J = 2
a = 1

H = let t = 1
	ax = CartesianIndices((-1:1, -1:1))
	hm = OffsetMatrix([zero(MMatrix{2, 2, typeof(t), 4}) for i in ax], ax)
	hm[1, 0][1, 2] = hm[0, 1][1, 2] = hm[-1, 0][2, 1] = hm[0, -1][2, 1] = hm[0, 0][1, 2] = hm[0, 0][2, 1] = t
	AutoBZCore.FourierSeries(SMatrix{2, 2, typeof(t), 4}.(hm), period = 1 / a)
end

tb_graphene = let t = 1.0
	ax = CartesianIndices((-2:2, -2:2))
	hm = OffsetMatrix([zero(MMatrix{2, 2, typeof(t), 4}) for i in ax], ax)
	hm[0, 0][1, 2] = hm[0, -1][1, 2] = hm[-1, 0][1, 2] = t
	hm[0, 0][2, 1] = hm[0, 1][2, 1] = hm[1, 0][2, 1] = t
	AutoBZCore.FourierSeries(SMatrix{2, 2, typeof(t), 4}.(hm), period = 1.0)
end
bz = load_bz(FBZ(), I(d))
function dos_graphene_exact(E::Real, t = oneunit(E))
	E = abs(E)
	x = abs(E / t)
	if x <= 1
		f = (1 + x)^2 - (x^2 - 1)^2 / 4
		2E / ((pi * t)^2 * sqrt(f)) * Elliptic.K(4x / f)
	elseif 1 < x < 3
		f = (1 + x)^2 - (x^2 - 1)^2 / 4
		2E / ((pi * t)^2 * sqrt(4x)) * Elliptic.K(f / 4x)
	else
		zero(inv(oneunit(t)))
	end
end

p0 = (; η = 1e-1, ω = 0.01)
greens_function(k, h_k, (; η, ω)) = -imag(tr(inv((ω + im * η) * I - h_k))) / (π * (2π)^d)
prototype = let k = AutoBZCore.FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end

integrand = FourierIntegralFunction(greens_function, H, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-3)


function dos_solver_iai(η, atol)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = atol), EvalCounter(IAI()))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		@SVector [getproperty(temp, :value), getproperty(getproperty(temp, :stats), :numevals)]
	end
end

function dos_solver_ptr(N, η)
	solver = init(AutoBZProblem(TrivialRep(), integrand, bz, p0), PTR(npt = N))
	ω -> begin
		solver.p = (; η, ω)
		temp = solve!(solver)

		getproperty(temp, :value)
	end
end


ηlist = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:3]))
tolerances = sort(unique([i * 10.0^(-j) for i in 1:10 for j in 1:5]))

