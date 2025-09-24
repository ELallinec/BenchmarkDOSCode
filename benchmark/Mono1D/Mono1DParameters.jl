using Pkg: Pkg
if !(Base.active_project() == "/workdir/ewen.lallinec/BCD/docs/Project.toml")
	Pkg.activate(".")           # reproducible environment included
end
push!(LOAD_PATH, "src")
using LinearAlgebra
using AutoBZCore
using FFTW
using Folds
using BenchmarkTools
using StaticArrays
using JLD2
#using DensityOfStates
using FourierSeriesEvaluators
#Common parameters
d = 1
#Hamiltonian and derivatives for BCD, LT and Smearing

C = MMatrix{1, 1, Float64, 1}[[1.0;;], [0.0;;], [1.0;;]]
H = FourierSeries(SMatrix{1, 1, Float64, 1}.(C), period = 1, offset = -2);
DH = AutoBZCore.HessianSeries(H);

function Ref(E, t)
	if abs(E) < 2t
		return 1 / (π * sqrt((2t)^2 - E^2))
	else
		return 0
	end
end


# IAI Parameters
HIAI = FourierSeries([1, 0, 1], period = 1, offset = -2);
tolerances = sort(unique([i * 10.0^j for j in -8:-1 for i in 1:0.1:10]), rev = true)
bz = load_bz(FBZ(d), I(d))

# IAI and Smearing parameters
ηlist = sort(unique([i * 10.0^j for j in -10:-1 for i in 1:0.1:10]), rev = true) #list of η for IAI and Smearing

# IAI and PTR solver 
p0 = (; η = 1e-1, ω = 0.01)
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = AutoBZCore.FourierSeriesEvaluators.period(H)
	greens_function(k, H(k), p0)
end

integrand = FourierIntegralFunction(greens_function, H, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-3)
prob = DOSProblem(H, float(zero(1.0)), bz)
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

