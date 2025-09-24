using Pkg: Pkg
Pkg.activate(".")           # reproducible environment included
Pkg.instantiate()           # install dependencies

using WannierIO

hrdat = read_w90_hrdat("data/svo_hr.dat")

Rmin, Rmax = extrema(hrdat.Rvectors)
Rsize = Tuple(Rmax .- Rmin .+ 1)
n, m = size(first(hrdat.H))

using StaticArrays, OffsetArrays

H_R = OffsetArray(
	Array{SMatrix{n, m, eltype(eltype(hrdat.H)), n * m}}(undef, Rsize...),
	map(:, Rmin, Rmax)...,
)
for (i, h, n) in zip(hrdat.Rvectors, hrdat.H, hrdat.Rdegens)
	H_R[CartesianIndex(Tuple(i))] = h / n
end

using FourierSeriesEvaluators, LinearAlgebra

h = FourierSeries(H_R, period = 1.0)

η = 1e-2                    # 10 meV (scattering amplitude)
ω_min = 10
ω_max = 15

using AutoBZCore
bz = load_bz(FBZ(), "data/svo.wout")
p0 = (; η, ω = (ω_min + ω_max) / 2) # initial parameters
greens_function(k, h_k, (; η, ω)) = tr(inv((ω + im * η) * I - h_k))
prototype = let k = FourierSeriesEvaluators.period(h)
	greens_function(k, h(k), p0)
end


integrand = FourierIntegralFunction(greens_function, h, prototype)
prob_dos = AutoBZProblem(TrivialRep(), integrand, bz, p0; abstol = 1e-3)

using HChebInterp

cheb_order = 15

function dos_solver(prob, alg)
	solver = init(prob, alg)
	ω -> begin
		solver.p = (; solver.p..., ω)
		solve!(solver).value
	end
end

dos_solver_ptr = dos_solver(prob_dos, PTR(; npt = 100))
@time greens_ptr = hchebinterp(dos_solver_ptr, ω_min, ω_max; atol = 1e-2, order = cheb_order)


dos_solver_iai = dos_solver(prob_dos, IAI())
@time greens_iai = hchebinterp(dos_solver_iai, ω_min, ω_max; atol = 1e-2, order = cheb_order)


#! BCD

prob = DOSProblem(h, float(zero(1.0)), bz)
function dos_solver_bcd(npt, η)
	solver = AutoBZCore.init(prob, AutoBZCore.BCD(; npt = npt, α = 0.1 / (2π), ΔE = 0.5, η = η))
	ω -> begin
		solver.domain = ω
		solve!(solver).value
	end
end
dosbcd = []

tabE = 11.5:0.1:14
@time dosiai = -imag.(dos_solver_iai.(tabE)) / π / det(bz.B)
@time dosiai_cheb = -imag.(greens_iai.(tabE)) / π / det(bz.B)
@time dosbcd = dos_solver_bcd(100, 1e-2).(tabE)
@time dosptr = -imag.(dos_solver_ptr.(tabE)) / π / det(bz.B)
@time dosptr_cheb = -imag.(greens_ptr.(tabE))/ π / det(bz.B)

using PyPlot
plot(tabE, dosptr, label = "PTR")
plot(tabE, dosptr_cheb, label = "PTR_Cheb")
plot(tabE, dosiai, label = "IAI")
plot(tabE, dosiai_cheb, label = "IAI_Cheb")
plot(tabE, dosbcd, label = "BCD")
legend()
