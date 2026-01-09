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
using QuadGK
a = 1
b = 3
H(k) = [a*k 0; 0 -b*k]


def(k, E, ΔE = 0.2) = -a*exp((-(a*k-E)^2)/(ΔE^2)) + b*exp((-(-b*k-E)^2)/(ΔE^2))
∂def(k, E, ΔE = 0.2)=ForwardDiff.derivative(k->def(k, E, ΔE), k)
bz(N, mul = 1) = (-N):(mul/N):N
function dos_bcd(Energies, N, α = 0.1, ΔE = 0.2, mul = 1)
	dos_values = zeros(ComplexF64, length(Energies))
	BZN=bz(N, mul)
	for (iE, E) in enumerate(Energies)

		for k in BZN
			G = inv(E*I - H(k + α*im*def(k, E, ΔE)))
			J = (1 + α*im*∂def(k, E, ΔE))
			dos_values[iE] += tr(G)*J
		end
	end
	return -imag.(dos_values)/(N*π)
end
ΔE=0.5
integrand(k, E, α = 0.1, ΔE = 0.2)=tr(inv((E+im*1e-16*(k==0.0))*I - H(k + α*im*def(k, E, ΔE)))*(1 + α*im*∂def(k, E, ΔE)))
Energies=range(-0.5, 0.5, 501)
@time temp=map(E->-imag(quadgk(k->integrand(k, E, 0.15, ΔE), -Inf, Inf, rtol = 1e-16)[1])/π, Energies);

N=100
mul = 100/17
@time tempbcd=dos_bcd(Energies, N, 0.2, ΔE, mul)
fig, ax = subplots()
ax.plot(Energies, 1/a .+ 1/b*ones(length(Energies)), lw = 5, color = :k, "-.", label = "Exact DOS")
ax.plot(Energies, temp, lw = 5, label = L"BCD \ N=\infty", color = :b)
ax.plot(Energies, tempbcd*mul, linestyle = :dashed, lw = 5, label = "BCD N=$(Int64(N/mul))", color = :r)
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xticks([-ΔE/4, 0, ΔE/4], [L"-\frac{ΔE}{4}", "0", L"\frac{ΔE}{4}"])
ax.set_yticks([1/b-1/a, 1/a+1/b], [L"\frac{1}{\beta}-\frac{1}{\gamma}", L"\frac{1}{\gamma}+\frac{1}{\beta}"])
ax.legend()
savefig("toymodeldosbcd.pdf")